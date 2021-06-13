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

    let (putative_phasing, contig_chunk_indices) = phase(&assembly, hic_mols, ccs, txg_barcodes, sex_contigs, &params);
    //phase(assembly, hic_mols, ccs, sex_contigs, &params);
    eprintln!("done phasing, writing phased vcf");
    output_phased_vcf(
        &kmers,
        &params,
        putative_phasing,
        &assembly,
        &contig_chunk_indices,
    );
    eprintln!("done");

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

fn phase(assembly: &Assembly, hic_mols: Mols, ccs_mols: Mols, txg_mols: Mols, sex_contigs: HashSet<i32>, params: &Params) -> (HashMap<i32, Vec<Option<bool>>>, HashMap<i32, Vec<(usize, usize)>>) {
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
    let mut contig_chunks: HashMap<i32, Vec<(usize, usize)>> = HashMap::new();
    let mut contig_phasing: HashMap<i32, Vec<Option<bool>>> = HashMap::new();


    for contig in 1..(assembly.contig_kmers.len()+1) {
        //if contig > 30 { break } // TODO remove
        if sex_contigs.contains(&(contig as i32)) { continue; }
        let length = *assembly.contig_sizes.get(&(contig as i32)).unwrap();
        if length <= params.min_contig_length { continue; }
        //let mut possible_positions: HashSet<usize> = HashSet::new();
        let kmer_positions = assembly.contig_kmers.get(&(contig as i32)).expect("please no");
        let mut putative_phasing: Vec<Option<bool>> = Vec::new();
        

        eprintln!("PHASING CONTIG {} id {} with {} kmer positions", assembly.contig_names[contig], contig, kmer_positions.len());

        let mut phase_blocks: HashMap<usize, (usize, usize)> = HashMap::new(); // map of phase block id to start, stop
        let mut current_phase_block_id: usize = 0;
        let mut position_phase_block : Vec<Option<usize>> = Vec::new();
        for _ in 0..kmer_positions.len() { // fill in phasings with unphased and we will come back and put in phasings as we... phase
            putative_phasing.push(None);
            position_phase_block.push(None);
        }

        
        let mut kmer_to_index: HashMap<i32, usize> = HashMap::new();
        for (index, (_pos, kmer)) in kmer_positions.iter().enumerate() {
            kmer_to_index.insert(Kmers::canonical_pair(*kmer), index);
        }

        
        let mut seeder: RandSeeder = RandSeeder::new(kmer_positions.len());
        let mut used_ccs_mols: BitSet = BitSet::new();
        let mut used_txg_mols: BitSet = BitSet::new();
        let mut deferred_seed: Option<usize> = None;
        let mut current_phase_block_start = 0;
        let mut current_phase_block_end = 0;
        let mut max_phase_block_id = 0;
        
        let mut kmer_phasing_consistency_counts: HashMap<i32, [u8; 4]> = HashMap::new();

        let mut no_counts_counter = 0;
        'outer_loop:
        loop { // loop over multiple phase blocks
            if let Some(seed_index) = deferred_seed { // going backwards if we have a deferred_seed
                eprintln!("continuing backwards at seed index {}", seed_index);
                deferred_seed = None;
                no_counts_counter = 0;
                let mut last_index = seed_index;
                for index in (0..seed_index).rev() { // going backwards
                    seeder.consume(index);
                    last_index = index;
                    let (position, kmer) = kmer_positions[index];
                    let canonical_kmer = Kmers::canonical_pair(kmer);
                    if let Some(counts) = kmer_phasing_consistency_counts.get(&canonical_kmer) {
                        let consistency = is_phasing_consistent(counts, &thresholds);
                        eprintln!("backwards kmer {}, position {}, index {}, counts {:?}, consistency {:?}", canonical_kmer, position, index, counts, consistency);
                        if consistency.is_consistent {
                            no_counts_counter = 0;
                            if let Some(overlapping_block) = position_phase_block[index] {
                                // make function to merge blocks and end
                                phase_blocks.insert(current_phase_block_id, (current_phase_block_start, current_phase_block_end));
                                let phasing = putative_phasing[index].unwrap();
                                // phasing, cis
                                // so !(phasing ^ cis) gives true for true/true and false/false and false otherwise
                                let cis = !(consistency.cis ^ phasing);
                                eprintln!("reverse merging blocks {} and {} in {} because overlapping kmer wants to be added in {} and has phase {} in its original block", current_phase_block_id, overlapping_block, cis, consistency.cis, phasing);
                                let (_new_block_id, _new_start, _new_end) = merge_phase_blocks(&mut phase_blocks, 
                                    &mut position_phase_block, &mut putative_phasing, 
                                    current_phase_block_id, overlapping_block, cis);
                                break;
                            }
                            current_phase_block_start = index;
                            position_phase_block[index] = Some(current_phase_block_id);
                            putative_phasing[index] = Some(consistency.cis);
                            add_kmer_and_update_phasing_consistency_counts(&mut kmer_phasing_consistency_counts, 
                                canonical_kmer, consistency.cis, &ccs_kmer_mols, &ccs_mols, &mut used_ccs_mols, &kmer_to_index, &kmer_positions, position, params.max_linked_read_dist);
                            add_kmer_and_update_phasing_consistency_counts(&mut kmer_phasing_consistency_counts, 
                                canonical_kmer, consistency.cis, &txg_kmer_mols, &txg_mols, &mut used_txg_mols, &kmer_to_index, &kmer_positions, position, params.max_linked_read_dist);
                        } else {
                            no_counts_counter += 1;
                        }
                    } else {
                        phase_blocks.insert(current_phase_block_id, (current_phase_block_start, current_phase_block_end));
                        eprintln!("backwards kmer {}, position {}, index {}, NOCOUNTS", canonical_kmer, position, index);
                        no_counts_counter += 1;
                        if no_counts_counter > 20 {
                            for backdex in ((index-20).max(0)..index).rev() {
                                let (_position, kmer) = kmer_positions[backdex];
                                let canonical_kmer = Kmers::canonical_pair(kmer);
                                if let Some(counts) = kmer_phasing_consistency_counts.get(&canonical_kmer) {
                                    eprintln!("\treaching backwards just to check position {}, index {} with {:?}", position, backdex, counts);
                                } else {
                                    eprintln!("\treaching backwards just to check position {}, index {} with NO COUNTS", position, backdex);
                                }
                            }
                            break;
                        }
                    }
                }
                phase_blocks.insert(current_phase_block_id, (current_phase_block_start, current_phase_block_end));
                eprintln!("backwards end index {}", current_phase_block_start);
                current_phase_block_id = max_phase_block_id + 1;
                max_phase_block_id += 1;
            } else {
                // going forwards, get new good seed
                let mut any = false;
                'seed_loop:
                    while let Some(seed_index) = seeder.next() {
                        any = true;
                        seeder.consume(seed_index);
                        let mut new_seed_bailout_count = 0;
                        let (position, kmer) = kmer_positions[seed_index]; // position is base position, index is the... index
                        let canonical_kmer = Kmers::canonical_pair(kmer.abs());
                        let kmer_consistency = *pairwise_kmer_consistency_counts.get(&canonical_kmer).unwrap_or(&0) as f32;
                        if !(kmer_consistency > min_seed_consistency && kmer_consistency < max_seed_consistency) {
                            seeder.consume(seed_index);
                            eprintln!("bad seed with {:?}", kmer_consistency);
                            continue;
                        } else {
                            no_counts_counter = 0;
                            new_seed_bailout_count = 0;
                            kmer_phasing_consistency_counts.clear(); // = HashMap::new();
                            eprintln!("found good seed {} with {:?} at seed index {}", canonical_kmer, kmer_consistency, seed_index);
                            deferred_seed = Some(seed_index.clone()); // start back here when done going forward
                            putative_phasing[seed_index] = Some(true);
                            add_kmer_and_update_phasing_consistency_counts(&mut kmer_phasing_consistency_counts, 
                                canonical_kmer, true, &ccs_kmer_mols, &ccs_mols, &mut used_ccs_mols, &kmer_to_index, &kmer_positions, position, params.max_linked_read_dist);
                            add_kmer_and_update_phasing_consistency_counts(&mut kmer_phasing_consistency_counts, 
                                canonical_kmer, true, &txg_kmer_mols, &txg_mols, &mut used_txg_mols, &kmer_to_index, &kmer_positions, position, params.max_linked_read_dist);
                            current_phase_block_start = seed_index;
                            current_phase_block_end = seed_index;
                            position_phase_block[seed_index] = Some(current_phase_block_id);

                            for index in (seed_index+1)..kmer_positions.len() { // going forward
                                seeder.consume(index);
                                let (position, kmer) = kmer_positions[index];
                                let canonical_kmer = Kmers::canonical_pair(kmer);
                                if let Some(counts) = kmer_phasing_consistency_counts.get(&canonical_kmer) {
                                    let consistency = is_phasing_consistent(counts, &thresholds);
                                    eprintln!("forwards kmer {}, position {}, index {}, counts {:?}, consistency {:?}", canonical_kmer, position, index, counts, consistency);
                                    if consistency.is_consistent {
                                        new_seed_bailout_count = 0;
                                        no_counts_counter = 0;
                                        if let Some(overlapping_block) = position_phase_block[index] {
                                            phase_blocks.insert(current_phase_block_id, (current_phase_block_start, current_phase_block_end));
                                            let phasing = putative_phasing[index].unwrap();
                                            // phasing, cis
                                             // so !(phasing ^ cis) gives true for true/true and false/false and false otherwise
                                            let cis = !(consistency.cis ^ phasing);
                                            eprintln!("forward merging blocks {} and {} in {} because overlapping kmer wants to be added in {} and has phase {} in its original block", 
                                                current_phase_block_id, overlapping_block, cis, consistency.cis, phasing);
                                            let (new_block_id, new_start, new_end) = merge_phase_blocks(&mut phase_blocks, 
                                                &mut position_phase_block, &mut putative_phasing, 
                                                current_phase_block_id, overlapping_block, cis);
                                            current_phase_block_id = new_block_id;
                                            assert!(new_start == current_phase_block_start);
                                            current_phase_block_start = new_start;
                                            current_phase_block_end = new_end;
                                            //write function to merge blocks and break 'seed_loop
                                            break 'seed_loop;
                                        } 
                                        putative_phasing[index] = Some(consistency.cis);
                                        add_kmer_and_update_phasing_consistency_counts(&mut kmer_phasing_consistency_counts, 
                                            canonical_kmer, consistency.cis, &ccs_kmer_mols, &ccs_mols, &mut used_ccs_mols, &kmer_to_index, &kmer_positions, position, params.max_linked_read_dist);
                                        add_kmer_and_update_phasing_consistency_counts(&mut kmer_phasing_consistency_counts, 
                                            canonical_kmer, consistency.cis, &txg_kmer_mols, &txg_mols, &mut used_txg_mols, &kmer_to_index, &kmer_positions, position, params.max_linked_read_dist);
                                        current_phase_block_end = index;
                                    } else {
                                        no_counts_counter += 1;
                                        new_seed_bailout_count += 1;
                                        eprintln!("checking bailout, new_seed_bailout_count {}, index {}, seed index {}", new_seed_bailout_count, index, seed_index);
                                        if new_seed_bailout_count > 10 && index - seed_index < 20 {
                                            eprintln!("FAILED SEED, do not pass go, do not collect 200$");
                                            for baddex in seed_index..(index + 1) {
                                                position_phase_block[baddex] = None;
                                                putative_phasing[baddex] = None;
                                                if baddex != seed_index {
                                                    seeder.unconsume(baddex);
                                                }
                                            }
                                            continue 'outer_loop;
                                        }
                                    }
                                } else {
                                    eprintln!("forward end kmer {}, position {}, index {}, NOCOUNTS", canonical_kmer, position, index);
                                    no_counts_counter += 1;
                                    new_seed_bailout_count += 1;
                                    if no_counts_counter > 20 {
                                        for fordex in (index+1).min(kmer_positions.len())..(index+20).min(kmer_positions.len()) {
                                            let (position, kmer) = kmer_positions[fordex];
                                            let canonical_kmer = Kmers::canonical_pair(kmer);
                                            if let Some(counts) = kmer_phasing_consistency_counts.get(&canonical_kmer) {
                                                eprintln!("\treaching forward just to check position {}, index {} with {:?}", position, fordex, counts);
                                            } else {
                                                eprintln!("\treaching forward just to check position {}, index {} with NO COUNTS", position, fordex);
                                            }
                                        }
                                        break 'seed_loop;
                                    }
                                }
                            }
                            eprintln!("end of contig");
                            break 'seed_loop;
                        }
                } // end forward seed loop
                // no more seeds
                if ! any {
                    let mut count_vec: Vec<(&usize, &(usize, usize))> = phase_blocks.iter().collect();
                    count_vec.sort_by(|a, b| b.1.cmp(a.1));
                    for (phase_block_id, (start, end)) in count_vec.iter() {
                        eprintln!("phase block {} goes from {}-{}", phase_block_id, kmer_positions[*start].0, kmer_positions[*end].0);
                    }
                    eprintln!("no more seeds, done with contig");
                    break 'outer_loop;
                }
                
            } // end forward/backward conditional 
        } // end phase block loop

        let phase_block_consistencies = get_phase_block_consistencies(&phase_blocks, &putative_phasing, &kmer_positions, &hic_mols, &hic_kmer_mols);
        
        let hic_thresholds = PhasingConsistencyThresholds{
            min_count: 20,
            min_percent: 0.75,
            minor_allele_fraction: 0.25,
        };
 


        let mut blocks: Vec<(&usize, &(usize, usize))> = phase_blocks.iter().collect();
        blocks.sort_by(|a, b| a.1.cmp(b.1));
        let mut new_blocks: Vec<(usize, usize)> = Vec::new();
        let (start, end) = blocks[0].1;
        let mut start = start;
        let mut end = end;
        for blockdex in 0..blocks.len() {
            let (_, new_end) = blocks[blockdex].1;
            end = new_end;
            if blockdex + 1 == blocks.len() {
                new_blocks.push((*start, *end));
            } else {
                let block_id1 = *blocks[blockdex].0;
                let block_id2 = *blocks[blockdex + 1].0;
                let counts = phase_block_consistencies.get(&(block_id1.min(block_id2), block_id1.max(block_id2))).unwrap_or(&[0;4]);
                let consistency = is_phasing_consistent(counts, &hic_thresholds);
                if !consistency.is_consistent {
                    new_blocks.push((*start, *end));
                    let (new_start, _) = blocks[blockdex + 1].1;
                    start = new_start;
                } else{
                    if !consistency.cis {
                        let (start, end) = *phase_blocks.get(&block_id2).expect("losing my mind if this fails");
                        for index in start..end {
                            if let Some(phase) = putative_phasing[index] {
                            putative_phasing[index] = Some(!phase);
                            }
                        }
                    }
                }
            }
        }

        eprintln!("final phase blocks for contig {}", contig);
        for (block_id, (start, end)) in new_blocks.iter().enumerate() {
            eprintln!("final phase block {} {}-{}", block_id, kmer_positions[*start].0, kmer_positions[*end].0);
        }
        contig_chunks.insert(contig as i32, new_blocks);
        contig_phasing.insert(contig as i32, putative_phasing);

    } // end contig loop

    let reader = fasta::Reader::from_file(Path::new(&params.assembly_fasta.to_string()))
        .expect("fasta not found");
    let mut writer = fasta::Writer::to_file(Path::new(&format!("{}/breaks.fa", params.output)))
        .expect("cannot open fasta writer");
    for record in reader.records() {
        let record = record.unwrap();
        let contig_name = record.id().to_string();
        let contig_id = assembly.contig_ids.get(&contig_name).unwrap();
        let kmer_positions = assembly.contig_kmers.get(contig_id).expect("don't even");


        let mut ranges: Vec<(usize, usize)> = Vec::new();
        if !contig_chunks.contains_key(contig_id) {
            eprintln!("contig has no chunks??? {}", contig_id);
            //let range = contig_chunks.entry(*contig_id).or_insert(Vec::new());
            ranges.push((0, *assembly.contig_sizes.get(contig_id).expect("bad things happen sometimes")));
        } else {
            let ranges_indices = contig_chunks.get(contig_id).unwrap();
            for (start, stop) in ranges_indices {
                ranges.push((kmer_positions[*start].0, kmer_positions[*stop].0));
            }
        }
        

        for (index, (start, stop)) in ranges.iter().enumerate() {
            let mut new_contig_name = contig_name.to_string();
            if ranges.len() > 0 {
                let list = vec![
                    new_contig_name,
                    (index + 1).to_string(),
                    start.to_string(),
                    stop.to_string(),
                ];
                new_contig_name = list.join("_");
            }
            let seq: TextSlice = &record.seq()[*start..*stop];
            let record = Record::with_attrs(&new_contig_name, None, &seq);
            writer
                .write_record(&record)
                .expect("could not write record");
        }
    } 
    (contig_phasing, contig_chunks)
}

fn allele(kmer: i32) -> Allele {
    match kmer.abs() % 2 == 0 {
        true => Allele::Ref,
        false => Allele::Alt,
    }
}


#[derive(Debug, Clone, Copy)]
enum Allele {
    Ref,
    Alt,
}

fn output_phased_vcf(
    kmers: &Kmers,
    params: &Params,
    contig_phasing: HashMap<i32, Vec<Option<bool>>>,
    assembly: &Assembly,
    contig_chunk_indices: &HashMap<i32, Vec<(usize, usize)>>,
) {
    let mut output = params.output.to_string();
    output.push_str("/phasing_breaks.vcf");
    let f = File::create(output).expect("Unable to create file");
    let mut f = BufWriter::new(f);
    let phasing: HashMap<i32, Vec<Option<bool>>> = contig_phasing;//HashMap::new();
    //for (contig, center) in best_centers.iter() {
    for contig in 1..assembly.contig_names.len() {
        let contig = &(contig as i32);

        
        //let contig_phasing = phasing.entry(*contig as i32).or_insert(Vec::new());
        let empty: Vec<Option<bool>> = Vec::new();
        let putative_phasing = phasing.get(contig).unwrap_or(&empty);
        
        let contig_name = &assembly.contig_names[*contig as usize];

        let chunks = contig_chunk_indices
            .get(contig)
            .expect("why do you hate me");
        let kmer_positions = assembly.contig_kmers.get(contig).expect("nooooo");
        let mut chunk_positions: Vec<(usize, usize)> = Vec::new();//contig_chunk_positions.get(contig).expect("noooo");
        for (start, end) in chunks.iter() {
            chunk_positions.push((kmer_positions[*start].0, kmer_positions[*end].0));
        }

        for (chunkdex, (left, right)) in chunks.iter().enumerate() {
            let left_pos = chunk_positions[chunkdex].0;
            let right_pos = chunk_positions[chunkdex].1;
            let chunk_name = vec![
                contig_name.to_string(),
                (chunkdex + 1).to_string(),
                left_pos.to_string(),
                right_pos.to_string(),
            ]
            .join("_");
            for ldex in *left..*right {
                //0..center.clusters[0].center.len() {
                // output is semi-vcf contig\tpos\t.\tREF\tALT\tqual\tfilter\tinfo\tformat\tsample

                let (pos, kmer) = kmer_positions[ldex];
                let reference;
                let alternate;
                let flip;
                match allele(kmer) {
                    Allele::Ref => {
                        reference = kmers.kmers.get(&kmer).unwrap().to_string();
                        alternate = kmers.kmers.get(&Kmers::pair(kmer)).unwrap().to_string();
                        flip = false;
                    }
                    Allele::Alt => {
                        reference = kmers.kmers.get(&Kmers::pair(kmer)).unwrap().to_string();
                        alternate = kmers.kmers.get(&kmer).unwrap().to_string();
                        flip = true;
                    }
                }


                let mut final_phasing: Vec<Option<bool>> = Vec::new();
                let genotype: String;
                if let Some(phase) = putative_phasing[ldex] {
                    if phase{
                        if !flip {
                            final_phasing.push(Some(true));
                            genotype = "0|1:60".to_string();
                        } else {
                            final_phasing.push(Some(false));
                            genotype = "1|0:60".to_string();
                        }
                    } else {
                        if !flip {
                            final_phasing.push(Some(false));
                            genotype = "1|0:60".to_string();
                        } else {
                            final_phasing.push(Some(true));
                            genotype = "0|1:60".to_string();
                        }
                    }
                } else {
                    final_phasing.push(None);
                    genotype = "./.:15".to_string();
                }

                
                let line_vec: Vec<String> = vec![
                    chunk_name.to_string(),
                    pos.to_string(),
                    ".".to_string(),
                    reference,
                    alternate,
                    ".".to_string(),
                    ".".to_string(),
                    ".".to_string(),
                    "GT:PQ".to_string(),
                    genotype.to_string(),
                ];
                let mut line = line_vec.join("\t");
                line.push_str("\n");
                f.write_all(line.as_bytes()).expect("Unable to write data");
            }
        }
    }
}

fn merge_phase_blocks(phase_blocks: &mut HashMap<usize, (usize, usize)>, 
    position_phase_block: &mut Vec<Option<usize>>, putative_phasing: &mut Vec<Option<bool>>, 
        phase_block1: usize, phase_block2: usize, cis: bool) -> (usize, usize, usize) {
    let new_start: usize;
    let new_end: usize;
    { // this scope is for the rust gods
        let (start1, end1) = phase_blocks.get(&phase_block1).unwrap();
        let (start2, end2) = phase_blocks.get(&phase_block2).unwrap();
        eprintln!("merging blocks {} and {} with positions {}-{} and {}-{}", phase_block1, phase_block2, start1, end1, start2, end2);
        new_start = *start1.min(start2);
        new_end = *end1.max(end2);
    }
    let canonical = phase_block1.min(phase_block2);
    let non_canonical = phase_block1.max(phase_block2);
    let (s1_remove, e1_remove) = phase_blocks.get(&non_canonical).unwrap();
    for index in *s1_remove..(e1_remove+1) {
        if let Some(phase) = putative_phasing[index] {
            position_phase_block[index] = Some(canonical);
            if !cis {
                putative_phasing[index] = Some(!phase);
            }
        }
    }
    eprintln!("removing block {}", non_canonical);
    { // why do i have to create random scopes?
        phase_blocks.remove(&non_canonical);
    }
    eprintln!("replacing with new block {} from {}-{}", canonical, new_start, new_end);
    phase_blocks.insert(canonical, (new_start, new_end));
    eprintln!("printing all blocks");
    for (id,(start, end)) in phase_blocks.iter() {
        eprintln!("\tblock {} from {}-{}", id, start, end);
    }

    (canonical, new_start, new_end)
}

fn get_phase_block_consistencies(phase_blocks: &HashMap<usize, (usize, usize)>, putative_phasing: &Vec<Option<bool>>, 
    kmer_positions: &Vec<(usize, i32)>, hic_mols: &Mols, hic_kmer_mols: &KmerMols) -> HashMap<(usize, usize), [u8; 4]> {
    let mut phase_block_consistencies: HashMap<(usize, usize), [u8; 4]> = HashMap::new();
    let mut blocks: Vec<usize> = Vec::new();
    
    
    let mut block_kmer_phasings: HashMap<usize, HashMap<i32, bool>> = HashMap::new();
    for (block_id, (start, end)) in phase_blocks.iter() {
        blocks.push(*block_id);
        let kmer_phasings = block_kmer_phasings.entry(*block_id).or_insert(HashMap::new());
        for index in *start..*end {
            let (_pos, kmer) = kmer_positions[index];
            let canonical_kmer = Kmers::canonical_pair(kmer);
            if let Some(phase) = putative_phasing[index] {
                kmer_phasings.insert(canonical_kmer, phase);
            }
        }
    }
    blocks.sort_by(|a, b| b.cmp(a));

   

    for phase_block1 in 0..blocks.len() {
        let phase_block1 = blocks[phase_block1];
        //let (start1, end1) = phase_blocks.get(&phase_block1).unwrap();

        let block1_phasing = block_kmer_phasings.get(&phase_block1).unwrap();
        for phase_block2 in (phase_block1 + 1)..phase_blocks.len() {
            let phase_block2 = blocks[phase_block2];
            //let (start2, end2) = phase_blocks[phase_block2];
            let block2_phasing = block_kmer_phasings.get(&phase_block2).unwrap();
            for (kmer1, phase1) in block1_phasing.iter() {
                let mut phase1 = *phase1;
                for mol in hic_kmer_mols.get_mols(*kmer1) {
                    for kmer2 in hic_mols.get_molecule_kmers(*mol) {
                        if let Some(phase2) = block2_phasing.get(kmer2) {
                            let mut phase2 = *phase2;
                            if kmer2 % 2 == 0 {
                                phase2 = !phase2;
                            }
                            let counts = phase_block_consistencies.entry((phase_block1, phase_block2)).or_insert([0;4]);
                            if phase1 && phase2 {
                                counts[0] += 1;
                            } else if !phase1 && !phase2 {
                                counts[1] += 1;
                            } else if phase1 && !phase2 {
                                counts[2] += 1;
                            } else {
                                counts[3] += 1;
                            }
                        }
                    }
                }
                for mol in hic_kmer_mols.get_mols(Kmers::pair(*kmer1)) {
                    let phase1 = !phase1;
                    for kmer2 in hic_mols.get_molecule_kmers(*mol) {
                        if let Some(phase2) = block2_phasing.get(kmer2) {
                            let mut phase2 = *phase2;
                            if kmer2 % 2 == 0 {
                                phase2 = !phase2;
                            }
                            let counts = phase_block_consistencies.entry((phase_block1, phase_block2)).or_insert([0;4]);
                            if phase1 && phase2 {
                                counts[0] += 1;
                            } else if !phase1 && !phase2 {
                                counts[1] += 1;
                            } else if phase1 && !phase2 {
                                counts[2] += 1;
                            } else {
                                counts[3] += 1;
                            }
                        }
                    }
                }
            }
        }
    }

    phase_block_consistencies
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
            
            let index = kmer_to_index.get(&canonical_kmer).unwrap_or(&0);
            
            let new_position = kmer_positions[*index].0 as i32;
            if (new_position - (current_position as i32)).abs() < max_distance as i32 {
                let mut counts = kmer_phasing_consistency_counts.entry(canonical_kmer).or_insert([0;4]);
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
            
            let index = kmer_to_index.get(&canonical_kmer).unwrap_or(&0);
            let new_position = kmer_positions[*index].0 as i32;
            if (new_position - (current_position as i32)).abs() < max_distance as i32 {
                let mut counts = kmer_phasing_consistency_counts.entry(canonical_kmer).or_insert([0;4]);
                increment_consistency_counts(!cis, *new_kmer, &mut counts);
            }
            
        }
        used.insert(*moldex);
    }

}

struct RandSeeder {
    shuffled_vec: Vec<usize>,
    used: BitSet,
    current_position: usize, // current_position in shuffled_vec
} 

impl RandSeeder {
    fn new(size: usize) -> RandSeeder {
        let seed: [u8; 32] = [4; 32]; // guarranteed to be random, chosen by fair dice roll https://xkcd.com/221/
        let mut rng: StdRng = SeedableRng::from_seed(seed);
        let mut shuffled_vec: Vec<usize> = (0..size).collect::<Vec<usize>>();
        shuffled_vec.shuffle(&mut rng);
        RandSeeder {
            used: BitSet::new(),
            shuffled_vec: shuffled_vec,
            current_position: 0,
        }
    }


    fn consume(&mut self, index: usize) {
        /*
        let shuffled_index = *self.reverse_index.get(&index).expect("please don't");
        let tmp = self.shuffled_vec[shuffled_index];
        if shuffled_index < self.vec_size {
            self.vec_size -= 1;   
        }
        self.shuffled_vec[shuffled_index] = self.shuffled_vec[self.vec_size];
        self.reverse_index.insert(self.shuffled_vec[shuffled_index], shuffled_index);
        self.shuffled_vec[self.vec_size] = tmp;
        self.reverse_index.insert(self.shuffled_vec[self.vec_size], self.vec_size);
        */
       
        self.used.insert(index);
        //eprintln!("consumed {} checking {}", index,self.used.contains(index));
    }

    fn unconsume(&mut self, index: usize) {
        self.used.remove(index);
        self.shuffled_vec.push(index);
    }

    fn next(&mut self) -> Option<usize> {
        
        while self.current_position < self.shuffled_vec.len() {
            let index = self.shuffled_vec[self.current_position];
            self.current_position += 1;
            if self.used.contains(index) {
                continue;
            }
            eprintln!("returning Some({}) and is it in used? {}", index, self.used.contains(index));
            return Some(index);
        }
        None
        /*
        if self.current_position >= self.vec_size {
            return None
        } else {
            Some(self.shuffled_vec[self.current_position - 1])
        }
        */
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

    for contig_id in 1..(assembly.contig_names.len()) {
        let size = assembly.contig_sizes.get(&(contig_id as i32)).expect("I am actually going crazy");
        let kmers = match assembly.contig_kmers.get(&(contig_id as i32)) {
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

    eprintln!("detecting sex contigs. mean kmer count is {} and mean paired kmer density is {}", avg_cov, avg_density);
    eprintln!("kmer_depth\thet_kmer_density\tcontig_id\tcontig_name\tcontig_length\tcontig_classification\tsex_contig_cov_cutoff\tsex_density_cutoff");
    for (depth, density, contig) in densities.iter() {
        let length = *assembly.contig_sizes.get(&(*contig as i32)).unwrap();
        let mut sex = "autosome";
        if length > params.min_contig_length && *depth < params.sex_contig_cov_cutoff * avg_cov 
            && *density < params.sex_contig_het_kmer_density_cutoff * avg_density {
                sex_contigs.insert(*contig as i32);
            sex = "sex";
        }
        eprintln!("{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}", depth, density, contig, 
            assembly.contig_names[*contig as usize], 
            length, sex,  params.sex_contig_cov_cutoff * avg_cov, 
            params.sex_contig_het_kmer_density_cutoff * avg_density);

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
    min_contig_length: usize,
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

    let sex_contig_het_kmer_density_cutoff = params.value_of("sex_contig_het_kmer_density_cutoff").unwrap_or("0.5");
    let sex_contig_het_kmer_density_cutoff = sex_contig_het_kmer_density_cutoff.to_string().parse::<f32>().unwrap();

    let sex_contig_cov_cutoff = params.value_of("sex_contig_cov_cutoff").unwrap_or("0.8");
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

    let min_contig_length = params.value_of("min_contig_length").unwrap_or("100000");
    let min_contig_length = min_contig_length.to_string().parse::<usize>().unwrap();

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
        min_contig_length: min_contig_length,
    }
}