#[macro_use]
extern crate clap;
extern crate bio;
extern crate disjoint_set;
extern crate hashbrown;
extern crate phasst_lib;
extern crate rand;
extern crate bit_set;

use bit_set::BitSet;

use bio::io::fasta;
use bio::io::fasta::Record;
use bio::utils::TextSlice;
use std::path::Path;

use phasst_lib::{
    load_assembly_kmers, load_hic, load_hifi, load_linked_read_barcodes, Assembly, Mols, Kmers, KmerMols,
};
use rayon::prelude::*;


use disjoint_set::DisjointSet;
use hashbrown::{HashMap, HashSet};

use clap::App;


use std::fs::File;
use std::io::{BufWriter, Write};
use std::io::BufReader;
use std::io::BufRead;
use std::io::Read;

use rand::Rng;
use rand::rngs::StdRng;
use rand::SeedableRng;
use rand::thread_rng;
use rand::seq::SliceRandom;


fn main() {
    println!("Welcome to phasst phase");
    let params = load_params();
    eprintln!("loading kmers");
    let kmers = Kmers::load_kmers(&params.het_kmers);
    //let (_variants, molecules) = load_molecule_kmers(&params.txg_mols, &params.hic_mols, &params.longread_mols, &kmers);
    eprintln!("loading hic kmers");
    let hic_mols = load_hic(Some(&params.hic_mols), &kmers, false);
    eprintln!("loading long reads");
    let ccs = load_hifi(Some(&params.ccs_mols), &kmers);
    eprintln!("loading linked reads");
    let txg_barcodes = load_linked_read_barcodes(Some(&params.txg_mols), &kmers);
    eprintln!("loading assembly kmers");
    let assembly = load_assembly_kmers(&params.assembly_kmers, &params.assembly_fasta, &kmers);

    let sex_contigs = detect_sex_contigs(&assembly, &params);
    phase(assembly, hic_mols, ccs, txg_barcodes, sex_contigs, &params);
    //phase(assembly, hic_mols, ccs, sex_contigs, &params);

}

fn get_mean_sd_pairwise_consistencies(pairwise_kmer_consistency_counts: &HashMap<i32, u32>) -> (f32, f32, f32, f32) {
    let mut sum = 0;
    for (_kmer, consistency) in pairwise_kmer_consistency_counts.iter() {
        sum += consistency;
    }
    let denom = pairwise_kmer_consistency_counts.len() as f32;
    let mean_phasing_consistency = (sum as f32)/denom;
    let mut sum_of_squared_diffs = 0.0;
    for (_kmer, consistency) in pairwise_kmer_consistency_counts.iter() {
        sum_of_squared_diffs += f32::powf(*consistency as f32 - mean_phasing_consistency, 2.0);
    }
    sum_of_squared_diffs /= denom;
    let sd_phasing_consistency = f32::sqrt(sum_of_squared_diffs);
    let min_seed_consistency = mean_phasing_consistency - sd_phasing_consistency;
    let max_seed_consistency = mean_phasing_consistency + sd_phasing_consistency;
    (mean_phasing_consistency, sd_phasing_consistency, min_seed_consistency, max_seed_consistency)
}

fn get_pairwise_consistencies(ccs_mols: &Mols) -> HashMap<(i32, i32), [u8; 4]> {
    let mut pairwise_consistencies: HashMap<(i32, i32), [u8; 4]> = HashMap::new();
    for ccs_mol in ccs_mols.get_molecules() {
        for k1dex in 0..ccs_mol.len() {
            let k1 = ccs_mol[k1dex].abs();
            for k2dex in k1dex..ccs_mol.len() {
                let k2 = ccs_mol[k2dex].abs();
                let mut key1 = k1.min(k2); // making the key in the hashtable canonical for which one is first and which one is second in the tuple
                let mut key2 = k1.max(k2);
                if key1 % 2 == 0 { // making the key in the hashtable canonical for which allele is used as the key
                    key1 = key1 - 1;
                }
                if key2 % 2 == 0 {
                    key2 = key2 - 1;
                }
                let counts = pairwise_consistencies.entry((key1, key2)).or_insert([0;4]);
                let k1_ref = k1 % 2 == 0; // which allele ref or alt, pairs are 1,2   3,4 etc
                let k2_ref = k2 % 2 == 0;
                if k1_ref && k2_ref {
                    counts[0] += 1;
                } else if !k1_ref && !k2_ref {
                    counts[1] += 1;
                } else if k1_ref && !k2_ref {
                    counts[2] += 1;
                } else {
                    counts[3] += 1;
                }
            }
        }
    }
    pairwise_consistencies
}

fn count_kmer_consistencies(pairwise_consistencies: &HashMap<(i32, i32), [u8;4]>, params: &Params) -> HashMap<i32, u32> {
    let mut pairwise_kmer_consistency_counts: HashMap<i32, u32> = HashMap::new();
    let thresholds = PhasingConsistencyThresholds{
        min_count: params.min_phasing_consistency_counts,
        min_percent: params.min_phasing_consistency_percent,
        minor_allele_fraction: params.min_minor_allele_fraction,
    };
    for ((k1, k2), counts) in pairwise_consistencies.iter() {
        let consistency = is_phasing_consistent(counts, &thresholds);
        {
            let k1_counts = pairwise_kmer_consistency_counts.entry(*k1).or_insert(0);
            if consistency.is_consistent {
                *k1_counts += 1;
            }
        }
        let k2_counts = pairwise_kmer_consistency_counts.entry(*k2).or_insert(0);
        if consistency.is_consistent {
            *k2_counts += 1;
        }
    }
    pairwise_kmer_consistency_counts
}

fn phase(assembly: Assembly, hic_mols: Mols, ccs_mols: Mols, txg_mols: Mols, sex_contigs: HashSet<i32>, params: &Params) {
//fn phase(assembly: Assembly, hic_mols: Mols, ccs_mols: Mols, sex_contigs: HashSet<i32>, params: &Params) {
    eprintln!("phasing");
    let hic_kmer_mols = hic_mols.get_kmer_mols();
    let ccs_kmer_mols = ccs_mols.get_kmer_mols();
    let txg_kmer_mols = txg_mols.get_kmer_mols();
    


    let pairwise_consistencies: HashMap<(i32, i32), [u8;4]> = get_pairwise_consistencies(&ccs_mols);

    // count kmer consistencies to make sure we seed on good kmers
    let pairwise_kmer_consistency_counts: HashMap<i32, u32> = count_kmer_consistencies(&pairwise_consistencies, &params);
    
    // get average consistency counts
    let (_mean_phasing_consistency, _sd_phasing_consistency, min_seed_consistency, max_seed_consistency) = 
        get_mean_sd_pairwise_consistencies(&pairwise_kmer_consistency_counts);
    
    let thresholds = PhasingConsistencyThresholds {
        min_count: params.min_phasing_consistency_counts,
        min_percent: params.min_phasing_consistency_percent,
        minor_allele_fraction: params.min_minor_allele_fraction,
    };
    

    for contig in 1..(assembly.contig_kmers.len()+1) {
        if contig > 1 { break } // TODO remove
        //let mut possible_positions: HashSet<usize> = HashSet::new();
        let kmer_positions = assembly.contig_kmers.get(&(contig as i32)).expect("please no");
        let mut putative_phasing: Vec<Option<bool>> = Vec::new();
        for _ in 0..kmer_positions.len() { // fill in phasings with unphased and we will come back and put in phasings as we... phase
            putative_phasing.push(None);
        }

        let mut phase_blocks: Vec<(usize, usize)> = Vec::new();
        
        let mut kmer_to_index: HashMap<i32, usize> = HashMap::new();
        for (index, (_pos, kmer)) in kmer_positions.iter().enumerate() {
            kmer_to_index.insert(Kmers::canonical_pair(*kmer), index);
        }







        
        let mut seeder: RandSeeder = RandSeeder::new(kmer_positions.len());
        let mut used_ccs_mols: BitSet = BitSet::new();
        let mut used_txg_mols: BitSet = BitSet::new();
        let mut deferred_seed: Option<usize> = None;

        'outer_loop:
            loop {
                let mut kmer_phasing_consistency_counts: HashMap<i32, [u8; 4]> = HashMap::new();
                if let Some(seed_index) = deferred_seed {
                    // going backwards TODO
                    eprintln!("continuing backwards at seed index {}", seed_index);
                    for index in (0..seed_index).rev() { // going backwards
                        let (position, kmer) = kmer_positions[index];
                        let canonical_kmer = Kmers::canonical_pair(kmer);
                        if let Some(counts) = kmer_phasing_consistency_counts.get(&canonical_kmer) {
                            let consistency = is_phasing_consistent(counts, &thresholds);
                            eprintln!("backwards kmer {}, index {}, counts {:?}, consistency {:?}", canonical_kmer, index, counts, consistency);
                            if consistency.is_consistent {
                                putative_phasing[index] = Some(consistency.cis);
                                add_kmer_and_update_phasing_consistency_counts(&mut kmer_phasing_consistency_counts, 
                                    canonical_kmer, consistency.cis, &ccs_kmer_mols, &ccs_mols, &mut used_ccs_mols, &kmer_to_index, &kmer_positions, position, params.max_linked_read_dist);
                                add_kmer_and_update_phasing_consistency_counts(&mut kmer_phasing_consistency_counts, 
                                    canonical_kmer, consistency.cis, &txg_kmer_mols, &txg_mols, &mut used_txg_mols, &kmer_to_index, &kmer_positions, position, params.max_linked_read_dist);
                                seeder.consume(index);
                            }
                        } else {
                            let current_phase_block = phase_blocks.len() - 1;
                            phase_blocks[current_phase_block].0 = index - 1;
                            eprintln!("backwards kmer {}, index {}, NOCOUNTS", canonical_kmer, index);
                            deferred_seed = None;
                            break;
                        }
                    }
                } else {
                    // going forwards, get new good seed
                    'seed_loop:
                        while let Some(seed_index) = seeder.next() {
                            let (position, kmer) = kmer_positions[seed_index]; // position is base position, index is the... index
                            let canonical_kmer = Kmers::canonical_pair(kmer.abs());
                            let kmer_consistency = *pairwise_kmer_consistency_counts.get(&canonical_kmer).unwrap_or(&0) as f32;
                            phase_blocks.push((seed_index, 0));
                            
                            if !(kmer_consistency > min_seed_consistency && kmer_consistency < max_seed_consistency) {
                                seeder.consume(seed_index);
                                eprintln!("bad seed with {:?}", kmer_consistency);
                                continue;
                            } else {
                                eprintln!("found good seed {} with {:?} at seed index {}", canonical_kmer, kmer_consistency, seed_index);
                                deferred_seed = Some(seed_index.clone()); // start back here when done going forward
                                putative_phasing[seed_index] = Some(true);
                                add_kmer_and_update_phasing_consistency_counts(&mut kmer_phasing_consistency_counts, 
                                    canonical_kmer, true, &ccs_kmer_mols, &ccs_mols, &mut used_ccs_mols, &kmer_to_index, &kmer_positions, position, params.max_linked_read_dist);
                                add_kmer_and_update_phasing_consistency_counts(&mut kmer_phasing_consistency_counts, 
                                    canonical_kmer, true, &txg_kmer_mols, &txg_mols, &mut used_txg_mols, &kmer_to_index, &kmer_positions, position, params.max_linked_read_dist);
                                for index in (seed_index+1)..kmer_positions.len() { // going forward
                                    let (position, kmer) = kmer_positions[index];
                                    let canonical_kmer = Kmers::canonical_pair(kmer);
                                    if let Some(counts) = kmer_phasing_consistency_counts.get(&canonical_kmer) {
                                        let consistency = is_phasing_consistent(counts, &thresholds);
                                        eprintln!("forwards kmer {}, index {}, counts {:?}, consistency {:?}", canonical_kmer, index, counts, consistency);

                                        if consistency.is_consistent {
                                            putative_phasing[index] = Some(consistency.cis);
                                            add_kmer_and_update_phasing_consistency_counts(&mut kmer_phasing_consistency_counts, 
                                                canonical_kmer, consistency.cis, &ccs_kmer_mols, &ccs_mols, &mut used_ccs_mols, &kmer_to_index, &kmer_positions, position, params.max_linked_read_dist);
                                            add_kmer_and_update_phasing_consistency_counts(&mut kmer_phasing_consistency_counts, 
                                                canonical_kmer, consistency.cis, &txg_kmer_mols, &txg_mols, &mut used_txg_mols, &kmer_to_index, &kmer_positions, position, params.max_linked_read_dist);
                                            
                                        }
                                        seeder.consume(index);

                                    } else {
                                        eprintln!("forwards kmer {}, index {}, NOCOUNTS", canonical_kmer, index);
                                        let current_phase_block = phase_blocks.len() - 1;
                                        phase_blocks[current_phase_block].1 = index - 1;
                                        break 'seed_loop;
                                    }
                                    
                                }
                                break;
                            }
                    }
                }
            }
  
    }

}

struct PhasingConsistencyThresholds {
    min_count: usize,
    min_percent: f32,
    minor_allele_fraction: f32,
}

fn increment_consistency_counts(phase: bool, allele: i32, counts: &mut [u8;4]) {
    let allele = allele.abs() % 2 == 1;
    if phase && allele {
        counts[0] += 1;
    } else if !phase && !allele { 
        counts[1] += 1;
    } else if phase && !allele {
        counts[2] += 1;
    } else {
        counts[3] += 1;
    }
}


fn add_kmer_and_update_phasing_consistency_counts(kmer_phasing_consistency_counts: &mut HashMap<i32, [u8;4]>, 
    kmer: i32, cis: bool, kmer_mols: &KmerMols, mols: &Mols, used: &mut BitSet, kmer_to_index: &HashMap<i32, usize>,
    kmer_positions: &Vec<(usize, i32)>, current_position: usize, max_distance: usize) {
    for moldex in kmer_mols.get_mols(kmer.abs()) { // loop over molecules which have kmer then loop over molecules that have the pair
        if used.contains(*moldex) {
            continue;
        }
        for new_kmer in mols.get_molecule_kmers(*moldex) {
            let canonical_kmer = Kmers::canonical_pair(*new_kmer);
            let mut counts = kmer_phasing_consistency_counts.entry(canonical_kmer).or_insert([0;4]);
            let index = kmer_to_index.get(&canonical_kmer).unwrap_or(&0);
            let new_position = kmer_positions[*index].0 as i32;
            if (new_position - current_position as i32).abs() < max_distance as i32 {
                increment_consistency_counts(cis, *new_kmer, &mut counts);
            }
            
        }
        used.insert(*moldex);
    }
    for moldex in kmer_mols.get_mols(Kmers::pair(kmer.abs())) {
        if used.contains(*moldex) {
            continue;
        }
        for new_kmer in mols.get_molecule_kmers(*moldex) {
            let canonical_kmer = Kmers::canonical_pair(*new_kmer);
            let mut counts = kmer_phasing_consistency_counts.entry(canonical_kmer).or_insert([0;4]);
            let index = kmer_to_index.get(&canonical_kmer).unwrap_or(&0);
            let new_position = kmer_positions[*index].0 as i32;
            if (new_position - current_position as i32).abs() < max_distance as i32 {
                increment_consistency_counts(!cis, *new_kmer, &mut counts);
            }
            
        }
        used.insert(*moldex);
    }

}

struct RandSeeder {
    seeded_rng: StdRng,
    shuffled_vec: Vec<usize>,
    vec_size: usize, // as we use up positions in the phasing process, we will move them to the end and shrink the vec_size
    current_position: usize, // current_position in shuffled_vec
    reverse_index: HashMap<usize, usize>, // map from position index to index in shuffled_vec so we can 
                                          // find them and move them to the back of shuffled_vec as they are used in the phasing process
} 

impl RandSeeder {
    fn new(size: usize) -> RandSeeder {
        let seed: [u8; 32] = [4; 32]; // guarranteed to be random, chosen by fair dice roll https://xkcd.com/221/
        let mut rng: StdRng = SeedableRng::from_seed(seed);
        let mut shuffled_vec: Vec<usize> = (0..size).collect::<Vec<usize>>();
        let mut reverse_index: HashMap<usize, usize> = HashMap::new();
        for (index, shuffled_index) in shuffled_vec.iter().enumerate() {
            reverse_index.insert(*shuffled_index, index);
        }
        shuffled_vec.shuffle(&mut rng);
        RandSeeder {
            seeded_rng: rng,
            shuffled_vec: shuffled_vec,
            vec_size: size,
            current_position: 0,
            reverse_index: reverse_index,
        }
    }


    fn consume(&mut self, index: usize) {
        let shuffled_index = *self.reverse_index.get(&index).expect("please don't");
        let tmp = self.shuffled_vec[shuffled_index];
        if shuffled_index < self.vec_size {
            self.vec_size -= 1;   
        }
        self.shuffled_vec[shuffled_index] = self.shuffled_vec[self.vec_size];
        self.reverse_index.insert(self.shuffled_vec[shuffled_index], shuffled_index);
        self.shuffled_vec[self.vec_size] = tmp;
        self.reverse_index.insert(self.shuffled_vec[self.vec_size], self.vec_size);
    } 

    fn next(&mut self) -> Option<usize> {
        self.current_position += 1;
        if self.current_position >= self.vec_size {
            return None
        } else {
            Some(self.shuffled_vec[self.current_position - 1])
        }
    }

    fn backtrack(&mut self) -> bool {
        if self.current_position > 0 {
            return false
        }
        self.current_position -= 1;
        let tmp = self.shuffled_vec[self.current_position];
        let pos = self.seeded_rng.gen_range(self.current_position + 1, self.vec_size);
        self.reverse_index.insert(tmp, pos);
        self.reverse_index.insert(self.shuffled_vec[pos], self.current_position);
        self.shuffled_vec[self.current_position] = self.shuffled_vec[pos];
        self.shuffled_vec[pos] = tmp;
        return true
    }
}

#[derive(Debug)]
struct PhasingConsistency {
    is_consistent: bool,
    cis: bool,
}

fn is_phasing_consistent(counts: &[u8;4], thresholds: &PhasingConsistencyThresholds) -> PhasingConsistency {
    let cis = (counts[0] + counts[1]) as f32;
    let trans = (counts[2] + counts[3]) as f32;
    let total = cis + trans;
    if cis > trans {
        let consistent_percentage = cis/total;
        if total > thresholds.min_count as f32 && consistent_percentage > thresholds.min_percent {
            let minor_allele = counts[0].min(counts[1]) as f32;
            if minor_allele/cis > thresholds.minor_allele_fraction {
                return PhasingConsistency{is_consistent: true, cis: true};
            }
        } 
    } else {
        let consistent_percentage = trans/total;
        if total > thresholds.min_count as f32 && consistent_percentage > thresholds.min_percent {
            let minor_allele = counts[2].min(counts[3]) as f32;
            if minor_allele/trans > thresholds.minor_allele_fraction {
                return PhasingConsistency{is_consistent: true, cis: false};
            }
        } 
    }
    PhasingConsistency{is_consistent: false, cis: false}
}



fn detect_sex_contigs(assembly: &Assembly, params: &Params) -> HashSet<i32> {
    let mut sex_contigs: HashSet<i32> = HashSet::new();
    let mut densities: Vec<(f32, f32, usize)> = Vec::new();
    let mut cov_sum: f32 = 0.0;
    let mut denom: f32 = 0.0;
    let mut density_sum: f32 = 0.0;


    //eprintln!("ok how many contigs are there in the assembly {}", assembly.molecules.len());

    for contig_id in 1..assembly.contig_names.len() {
        let size = assembly.contig_sizes.get(&(contig_id as i32)).expect("I am actually going crazy");
        let kmers = match assembly.molecules.get(&(contig_id as i32)) {
            Some(x) => x.len(),
            None => 0,
        };
        //eprintln!("contig_id {}",contig_id);
        densities.push((params.contig_kmer_cov[contig_id], (kmers as f32)/(*size as f32), contig_id));
        let size = *size as f32;
        cov_sum += params.contig_kmer_cov[contig_id] * size;
        denom += size;
        density_sum += kmers as f32;
    }

    let avg_cov = cov_sum / denom;
    let avg_density = density_sum / denom;

    //eprintln!("kmer_depth\thet_kmer_density\tcontig_id\tcontig_name\tcontig_length\tcontig_classification\tsex_contig_cov_cutoff\tsex_density_cutoff");
    for (depth, density, contig) in densities.iter() {
        
        if *depth < params.sex_contig_cov_cutoff * avg_cov 
            && *density < params.sex_contig_het_kmer_density_cutoff * avg_density {
                sex_contigs.insert(*contig as i32);
            
        }
        //eprintln!("{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}", depth, density, contig, 
        //    assembly.contig_names[*contig as usize], 
        //    assembly.contig_sizes.get(&(*contig as i32)).unwrap(), 
        //    sex,  params.sex_contig_cov_cutoff * avg_cov, 
        //    params.sex_contig_het_kmer_density_cutoff * avg_density);

    }
    sex_contigs
}


#[derive(Clone)]
struct Params {
    het_kmers: String,
    txg_mols: Vec<String>,
    hic_mols: Vec<String>,
    ccs_mols: Vec<String>,
    contig_kmer_cov: Vec<f32>,
    output: String,
    assembly_kmers: String,
    assembly_fasta: String,
    threads: usize,
    seed: u8,
    ploidy: usize,
    sex_contig_het_kmer_density_cutoff: f32,
    sex_contig_cov_cutoff: f32,
    restarts: u32,
    min_hic_links: u32,
    min_minor_allele_fraction: f32,
    min_phasing_consistency_counts: usize,
    min_phasing_consistency_percent: f32,
    break_window: usize,
    max_linked_read_dist: usize,
}

fn load_params() -> Params {
    let yaml = load_yaml!("params.yml");
    let params = App::from_yaml(yaml).get_matches();

    let het_kmers = params.value_of("het_kmers").unwrap();
    let output = params.value_of("output").unwrap();


    let mut txg_mols: Vec<String> = Vec::new();
    match params.value_of("linked_read_mols") {
        Some(txg_fofn) => {
            let f = File::open(txg_fofn).expect("Unable to open txg fofn");
            let f = BufReader::new(f);

            for line in f.lines() {
                let line = line.expect("Unable to read txg fofn line");
                txg_mols.push(line.to_string());
            }
        },
        None => (),
    }


    let mut hic_mols: Vec<String> = Vec::new();
    match params.value_of("hic_mols") {
        Some(hic_fofn) => {
            let f = File::open(hic_fofn).expect("Unable to open hic fofn");
            let f = BufReader::new(f);

            for line in f.lines() {
                let line = line.expect("Unable to read txg fofn line");
                hic_mols.push(line.to_string());
            }
        },
        None => (),
    }
    

    let mut ccs_mols: Vec<String> = Vec::new();
    match params.value_of("long_read_mols") {
        Some(ccs_fofn) => {
            let f = File::open(ccs_fofn).expect("Unable to open hic fofn");
            let f = BufReader::new(f);

            for line in f.lines() {
                let line = line.expect("Unable to read txg fofn line");
                ccs_mols.push(line.to_string());
            }
        },
        None => (),
    }


    let mut contig_kmer_cov: Vec<f32> = Vec::new();
    contig_kmer_cov.push(0.0); // contig ids are 1 indexed
    let contig_kmer_cov_file = params.value_of("contig_kmer_depths").unwrap();
    let f = File::open(contig_kmer_cov_file).expect("Unable to open contig_kmer_depths file");
    let f = BufReader::new(f);

    for line in f.lines() {
        let line = line.expect("unable to read contig_kmer_depths line");
        let toks: Vec<&str> = line.split("\t").collect();
        contig_kmer_cov.push(toks[1].to_string().parse::<f32>().expect("cannot parse float in contig_kmer_depths file"));
    }

    let threads = params.value_of("threads").unwrap_or("1");
    let threads = threads.to_string().parse::<usize>().unwrap();

    let seed = params.value_of("seed").unwrap_or("4"); // 4 is guarranteed random by dice roll https://xkcd.com/221/
    let seed = seed.to_string().parse::<u8>().unwrap();

    let restarts = params.value_of("restarts").unwrap_or("10");
    let restarts = restarts.to_string().parse::<u32>().unwrap();

    let assembly_kmers = params.value_of("assembly_kmers").unwrap();
    let assembly_fasta = params.value_of("assembly_fasta").unwrap();

    let ploidy = params.value_of("ploidy").unwrap_or("2");
    let ploidy = ploidy.to_string().parse::<usize>().unwrap();

    let sex_contig_het_kmer_density_cutoff = params.value_of("sex_contig_het_kmer_density_cutoff").unwrap_or("0.1");
    let sex_contig_het_kmer_density_cutoff = sex_contig_het_kmer_density_cutoff.to_string().parse::<f32>().unwrap();

    let sex_contig_cov_cutoff = params.value_of("sex_contig_cov_cutoff").unwrap_or("0.65");
    let sex_contig_cov_cutoff = sex_contig_cov_cutoff.to_string().parse::<f32>().unwrap();
        
    let min_hic_links = params.value_of("min_hic_links").unwrap_or("4");
    let min_hic_links = min_hic_links.to_string().parse::<u32>().unwrap();

    let break_window = params.value_of("break_window").unwrap_or("500");
    let break_window = break_window.to_string().parse::<usize>().unwrap();
    eprintln!("break window {}", break_window);

    let min_minor_allele_fraction = params.value_of("min_minor_allele_fraction").unwrap_or("0.15");
    let min_minor_allele_fraction = min_minor_allele_fraction.to_string().parse::<f32>().unwrap();

    let min_phasing_consistency_counts = params.value_of("min_phasing_consistency_counts").unwrap_or("8");
    let min_phasing_consistency_counts = min_phasing_consistency_counts.to_string().parse::<usize>().unwrap(); 

    let min_phasing_consistency_percent = params.value_of("min_phasing_consistency_percent").unwrap_or("0.9");
    let min_phasing_consistency_percent = min_phasing_consistency_percent.to_string().parse::<f32>().unwrap();

    let max_linked_read_dist = params.value_of("min_linked_read_dist").unwrap_or("150000");
    let max_linked_read_dist = max_linked_read_dist.to_string().parse::<usize>().unwrap();

    Params {
        het_kmers: het_kmers.to_string(),
        output: output.to_string(),
        txg_mols: txg_mols,
        hic_mols: hic_mols,
        ccs_mols: ccs_mols,
        contig_kmer_cov: contig_kmer_cov,
        assembly_kmers: assembly_kmers.to_string(),
        assembly_fasta: assembly_fasta.to_string(),
        threads: threads,
        seed: seed,
        sex_contig_het_kmer_density_cutoff: sex_contig_het_kmer_density_cutoff,
        sex_contig_cov_cutoff: sex_contig_cov_cutoff,
        restarts: restarts,
        min_minor_allele_fraction: min_minor_allele_fraction,
        min_phasing_consistency_counts: min_phasing_consistency_counts,
        min_phasing_consistency_percent: min_phasing_consistency_percent,
        ploidy: ploidy,
        min_hic_links: min_hic_links,
        break_window: break_window,
        max_linked_read_dist: max_linked_read_dist,
    }
}