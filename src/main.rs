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
    //let hic_mols = load_hic(Some(&params.hic_mols), &kmers, false); TODO ADD BACK
    eprintln!("loading long reads");
    let ccs = load_hifi(Some(&params.ccs_mols), &kmers);
    eprintln!("loading linked reads");
    // let txg_barcodes = load_linked_read_barcodes(Some(&params.txg_mols), &kmers); TODO ADD BACK
    eprintln!("loading assembly kmers");
    let assembly = load_assembly_kmers(&params.assembly_kmers, &params.assembly_fasta, &kmers);
    eprintln!("detecting sex contigs");
    let sex_contigs = detect_sex_contigs(&assembly, &ccs, &params, &kmers);
    /*
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
    */
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

fn get_pairwise_consistencies(ccs_mols: &Mols, assembly: &Assembly, any_number: bool) -> HashMap<(i32, i32), [u32; 4]> {
    let mut pairwise_consistencies: HashMap<(i32, i32), [u32; 4]> = HashMap::new();
    for (mol_id, ccs_mol) in ccs_mols.get_molecules().enumerate() {
        for k1dex in 0..ccs_mol.len() {
            let k1 = ccs_mol[k1dex].abs();
            if let Some((contig1, pos1)) = kmer_contig_position(k1, assembly, any_number){
                for k2dex in (k1dex+1)..ccs_mol.len() {
                    let k2 = ccs_mol[k2dex].abs();
                    if let Some((contig2, pos2)) = kmer_contig_position(k2, assembly, any_number) {
                        if contig1 != contig2 || pos1.max(pos2) - pos1.min(pos2) > 50000 {
                            continue;
                        }
                        let key1 = Kmers::canonical_pair(k1.abs().min(k2.abs()));
                        let key2 = Kmers::canonical_pair(k1.abs().max(k2.abs()));

                        let counts = pairwise_consistencies.entry((key1, key2)).or_insert([0;4]);

                        let k1_ref = k1.abs().min(k2.abs()) % 2 == 0; // which allele ref or alt, pairs are 1,2   3,4 etc
                        let k2_ref = k1.abs().max(k2.abs()) % 2 == 0;
                        if k1_ref && k2_ref {
                            counts[0] += 1;
                        } else if !k1_ref && !k2_ref {
                            counts[1] += 1;
                        } else if k1_ref && !k2_ref {
                            counts[2] += 1;
                        } else {
                            counts[3] += 1;
                        }
                        if ( pos1 == 15769  || pos1 == 19719 ||  pos1 == 20236 ) && (pos2 == 15769  || pos2 == 19719 ||  pos2 == 20236) {
                            eprintln!("lets dig in. pos1-pos2 {}-{} k1_ref-k2_ref {}-{} mol_id {} current counts {:?}", pos1, pos2, k1_ref, k2_ref, mol_id, counts);
                        }
                    }
                    
                    
                    
                }
            }
            
        }
    }
    pairwise_consistencies
}

fn kmer_contig_position(kmer: i32, assembly: &Assembly, any: bool) -> Option<(i32, usize)> {
    if let Some((contig_id, number_seen, _order, position)) = assembly.variants.get(&kmer.abs()) {
        if any || *number_seen == 1 {
            return Some((*contig_id, *position));
        }
    } else if let Some((contig_id, number_seen, _order, position)) = assembly.variants.get(&Kmers::pair(kmer.abs())) {
        if any || *number_seen == 1 {
            return Some((*contig_id, *position));
        }
    }
    None
}

fn count_kmer_consistencies(pairwise_consistencies: &HashMap<(i32, i32), [u32;4]>, params: &Params) -> HashMap<i32, u32> {
    let mut pairwise_kmer_consistency_counts: HashMap<i32, u32> = HashMap::new();
    let thresholds = PhasingConsistencyThresholds{
        min_count: params.min_phasing_consistency_counts,
        min_percent: params.min_phasing_consistency_percent,
        minor_allele_fraction: params.min_minor_allele_fraction,
    };
    
    for ((k1, k2), counts) in pairwise_consistencies.iter() {
        let consistency = is_phasing_consistent(counts, &thresholds, false);
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
    


    let pairwise_consistencies: HashMap<(i32, i32), [u32;4]> = get_pairwise_consistencies(&ccs_mols, assembly, false);

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
        //if contig > 1 { break } // TODO remove
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
        
        let mut kmer_phasing_consistency_counts: HashMap<i32, [u32; 4]> = HashMap::new();

        let mut no_counts_counter = 0;
        'outer_loop:
        loop { // loop over multiple phase blocks
            if let Some(seed_index) = deferred_seed { // going backwards if we have a deferred_seed
                //eprintln!("continuing backwards in phase block {} at seed index {}", current_phase_block_id, seed_index);
                deferred_seed = None;
                no_counts_counter = 0;
                let mut last_index = seed_index;
                for index in (0..seed_index).rev() { // going backwards
                    seeder.consume(index);
                    last_index = index;
                    let (position, kmer) = kmer_positions[index];
                    let canonical_kmer = Kmers::canonical_pair(kmer);
                    if let Some(counts) = kmer_phasing_consistency_counts.get(&canonical_kmer) {
                        let consistency = is_phasing_consistent(counts, &thresholds, false);
                        //eprintln!("backwards kmer {}, position {}, index {}, counts {:?}, consistency {:?}", canonical_kmer, position, index, counts, consistency);
                        if consistency.is_consistent {
                            no_counts_counter = 0;
                            if let Some(overlapping_block) = position_phase_block[index] {
                                // make function to merge blocks and end
                                phase_blocks.insert(current_phase_block_id, (current_phase_block_start, current_phase_block_end));
                                let phasing = putative_phasing[index].unwrap();
                                // phasing, cis
                                // so !(phasing ^ cis) gives true for true/true and false/false and false otherwise
                                let cis = !(consistency.cis ^ phasing);
                                //eprintln!("reverse merging blocks {} and {} in {} because overlapping kmer wants to be added in {} and has phase {} in its original block", current_phase_block_id, overlapping_block, cis, consistency.cis, phasing);
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
                        //phase_blocks.insert(current_phase_block_id, (current_phase_block_start, current_phase_block_end));
                        //eprintln!("backwards kmer {}, position {}, index {}, NOCOUNTS", canonical_kmer, position, index);
                        no_counts_counter += 1;
                        if no_counts_counter > 20 {
                            for backdex in ((index-20).max(0)..index).rev() {
                                let (_position, kmer) = kmer_positions[backdex];
                                let canonical_kmer = Kmers::canonical_pair(kmer);
                                /*
                                if let Some(counts) = kmer_phasing_consistency_counts.get(&canonical_kmer) {
                                    eprintln!("\treaching backwards just to check position {}, index {} with {:?}", position, backdex, counts);
                                } else {
                                    eprintln!("\treaching backwards just to check position {}, index {} with NO COUNTS", position, backdex);
                                }
                                */
                            }
                            break;
                        }
                    }
                }
                let total = (current_phase_block_end-current_phase_block_start) as f32;
                let mut phased = 0.0;
                for index in current_phase_block_start..current_phase_block_end {
                    if let Some(_) = putative_phasing[index] { phased += 1.0; }
                }
                if phased > 200.0 {
                    phase_blocks.insert(current_phase_block_id, (current_phase_block_start, current_phase_block_end));
                    eprintln!("backwards end, phase block {} goes from {}-{} indices which is {}-{} bases, length {}, {}% phased", current_phase_block_id, 
                        current_phase_block_start, current_phase_block_end, 
                        kmer_positions[current_phase_block_start].0, kmer_positions[current_phase_block_end].0,
                        kmer_positions[current_phase_block_end].0 - kmer_positions[current_phase_block_start].0, phased/total);
                } else {
                    eprintln!("backwards end, but this was a sucky phase block and we are dropping it, id {}, from {}-{} indices which is {}-{} length {} and had {}% phased", current_phase_block_id, 
                        current_phase_block_start, current_phase_block_end, 
                        kmer_positions[current_phase_block_start].0, kmer_positions[current_phase_block_end].0,
                        kmer_positions[current_phase_block_end].0 - kmer_positions[current_phase_block_start].0, phased/total);
                }
                
                
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
                            eprintln!("starting phase block {} found good seed {} with {:?} at seed index {}", current_phase_block_id, canonical_kmer, kmer_consistency, seed_index);
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
                                    let consistency = is_phasing_consistent(counts, &thresholds, false);
                                    //eprintln!("forwards kmer {}, position {}, index {}, counts {:?}, consistency {:?}", canonical_kmer, position, index, counts, consistency);
                                    if consistency.is_consistent {
                                        new_seed_bailout_count = 0;
                                        no_counts_counter = 0;
                                        if let Some(overlapping_block) = position_phase_block[index] {
                                            phase_blocks.insert(current_phase_block_id, (current_phase_block_start, current_phase_block_end));
                                            let phasing = putative_phasing[index].unwrap();
                                            // phasing, cis
                                             // so !(phasing ^ cis) gives true for true/true and false/false and false otherwise
                                            let cis = !(consistency.cis ^ phasing);
                                            //eprintln!("forward merging blocks {} and {} in {} because overlapping kmer wants to be added in {} and has phase {} in its original block", 
                                            //    current_phase_block_id, overlapping_block, cis, consistency.cis, phasing);
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
                                        //eprintln!("checking bailout, new_seed_bailout_count {}, index {}, seed index {}", new_seed_bailout_count, index, seed_index);
                                        if new_seed_bailout_count > 10 && index - seed_index < 20 {
                                            eprintln!("FAILED SEED, do not pass go, do not collect 200$");
                                            deferred_seed = None;
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
                                    //eprintln!("forward end kmer {}, position {}, index {}, NOCOUNTS", canonical_kmer, position, index);
                                    no_counts_counter += 1;
                                    new_seed_bailout_count += 1;
                                    if new_seed_bailout_count > 10 && index - seed_index < 20 {
                                        eprintln!("FAILED SEED, do not pass go, do not collect 200$");
                                        deferred_seed = None;
                                        for baddex in seed_index..(index + 1) {
                                            position_phase_block[baddex] = None;
                                            putative_phasing[baddex] = None;
                                            if baddex != seed_index {
                                                seeder.unconsume(baddex);
                                            }
                                        }
                                        continue 'outer_loop;
                                    }
                                    if no_counts_counter > 20 {
                                        for fordex in (index+1).min(kmer_positions.len())..(index+20).min(kmer_positions.len()) {
                                            let (position, kmer) = kmer_positions[fordex];
                                            let canonical_kmer = Kmers::canonical_pair(kmer);
                                            /*
                                            if let Some(counts) = kmer_phasing_consistency_counts.get(&canonical_kmer) {
                                                eprintln!("\treaching forward just to check position {}, index {} with {:?}", position, fordex, counts);
                                            } else {
                                                eprintln!("\treaching forward just to check position {}, index {} with NO COUNTS", position, fordex);
                                            }
                                            */
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
                    //for (phase_block_id, (start, end)) in count_vec.iter() {
                    //    eprintln!("phase block {} goes from {}-{}", phase_block_id, kmer_positions[*start].0, kmer_positions[*end].0);
                    //}
                    eprintln!("no more seeds, done with contig");
                    break 'outer_loop;
                }
                
            } // end forward/backward conditional 
        } // end phase block loop

        let mut phase_block_indices: HashMap<usize, Vec<usize>> = HashMap::new();
        for (phase_block_id, (start, end)) in phase_blocks.iter() {
            let indices = phase_block_indices.entry(*phase_block_id).or_insert(Vec::new());
            for i in *start..(*end + 1) {
                indices.push(i);
            }
        }

        let phase_block_consistencies = get_phase_block_consistencies(&phase_block_indices, &putative_phasing, &kmer_positions, &hic_mols, &hic_kmer_mols);
        
        //eprintln!("phase_block_consistencies.len() {}", phase_block_consistencies.len());
        let hic_thresholds = PhasingConsistencyThresholds{
            min_count: 20,
            min_percent: 0.75,
            minor_allele_fraction: 0.25,
        };

        let mut ordered_consistencies: Vec<((usize, usize), [u32; 4])> = Vec::new();
        for ((block_id1, block_id2), counts) in phase_block_consistencies.iter() {
            let mut copycounts: [u32;4] = [0;4];
            for i in 0..4 { copycounts[i] = counts[i]; }
            ordered_consistencies.push(((*block_id1, *block_id2), copycounts));
        }
        ordered_consistencies.sort_by(|a, b| b.1.iter().sum::<u32>().cmp(&a.1.iter().sum::<u32>()));



    
        // need to decide if we want to merge then recompare or just have all vs all then make connected component of phasing consistent stuff
        let mut any_merged = true; // TODO change to true, just turning this loop off for now.
        
        while any_merged {
            any_merged = false;
            for ((block_id1, block_id2), counts) in ordered_consistencies.iter() {
                let consistency = is_phasing_consistent(counts, &hic_thresholds, false);
                if consistency.is_consistent {
                    // do some merge process involving phase_blocks, putative_phasing 
                    eprintln!("merging blocks {} and {} with counts {:?} with sizes {} and {}", 
                        block_id1, block_id2, counts, phase_block_indices.get(block_id1).expect("blame richard").len(), 
                        phase_block_indices.get(block_id2).expect("blame richard").len());
                    if !consistency.cis {
                        for index in phase_block_indices.get(block_id2).expect("i expected otherwise") {
                            if let Some(phase) = putative_phasing[*index] { putative_phasing[*index] = Some(!phase); }
                        }
                    }
                    let mut hodler: Vec<usize> = Vec::new();
                    for i in phase_block_indices.get(block_id2).expect("blame sangjin") {
                        hodler.push(*i);
                    }
                    let indices = phase_block_indices.get_mut(block_id1).expect("nope");
                    for i in hodler {
                        indices.push(i);
                    }
                    phase_block_indices.remove(block_id2);
                    eprintln!("resulting block {} with length {}", block_id1, phase_block_indices.get(block_id1).expect("blame richard3").len());
                    
                    // recalculate phase_block_consistencies
                    let phase_block_consistencies = get_phase_block_consistencies(&phase_block_indices, &putative_phasing, &kmer_positions, &hic_mols, &hic_kmer_mols);
                    ordered_consistencies.clear();
                    for ((block_id1, block_id2), counts) in phase_block_consistencies.iter() {
                        let mut copycounts: [u32;4] = [0;4];
                        for i in 0..4 { copycounts[i] = counts[i]; }
                        ordered_consistencies.push(((*block_id1, *block_id2), copycounts));
                    }
                    ordered_consistencies.sort_by(|a, b| b.1.iter().sum::<u32>().cmp(&a.1.iter().sum::<u32>()));
                    any_merged = true;
                    break;
                } else {
                    eprintln!("NOT merging blocks {} and {} with counts {:?} with sizes {} and {}", 
                        block_id1, block_id2, counts, phase_block_indices.get(block_id1).expect("blame richard").len(), 
                            phase_block_indices.get(block_id2).expect("blame richard").len());
                    let mut min1 = usize::MAX;
                    let mut min2 = usize::MAX;
                    let mut max1 = 0;
                    let mut max2 = 0;
                    for index in phase_block_indices.get(block_id1).unwrap().iter() {

                        min1 = min1.min(kmer_positions[*index].0);
                        max1 = max1.max(kmer_positions[*index].0);
                    }
                    for index in phase_block_indices.get(block_id2).unwrap().iter() {

                        min2 = min2.min(kmer_positions[*index].0);
                        max2 = max2.max(kmer_positions[*index].0);
                    }
                    eprintln!("\tout of curiosity the position ranges for block {} and {} are {}-{} and {}-{}", block_id1, block_id2, min1, max1, min2, max2);
                    
                }

            }
        }

        for (block_id, block_indices) in phase_block_indices.iter() {
            eprintln!("contig {} final block id {} with size {}", contig, block_id, block_indices.len());
        }
        

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
                let start_pos = kmer_positions[*start].0 - kmer_positions[*start].0.min(30);
                let length = assembly.contig_sizes.get(contig_id).unwrap();
                let end_pos = kmer_positions[*stop].0 + (length - kmer_positions[*stop].0).min(30);
                ranges.push((start_pos, end_pos));
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
        true => Allele::Alt,
        false => Allele::Ref,
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

        let kmer_positions = assembly.contig_kmers.get(contig).expect("nooooo");
        
        //let contig_phasing = phasing.entry(*contig as i32).or_insert(Vec::new());
        let mut empty: Vec<Option<bool>> = Vec::new();
        if !phasing.contains_key(contig) {
            for _ in kmer_positions.iter() {
                empty.push(None);
            }   
        }
        let putative_phasing = phasing.get(contig).unwrap_or(&empty);
       
        let contig_name = &assembly.contig_names[*contig as usize];
        let mut chunk_positions: Vec<(usize, usize)> = Vec::new();
        let mut chunks: Vec<(usize, usize)> = Vec::new();
        if !contig_chunk_indices.contains_key(contig) {
            chunk_positions.push((0, *assembly.contig_sizes.get(contig).expect("really?")));
            chunks.push((0, kmer_positions.len()));
            
        } else {
            let chunk_indices = contig_chunk_indices
            .get(contig)
            .expect("why do you hate me");
            for (start, end) in chunk_indices.iter() {
                chunk_positions.push((kmer_positions[*start].0, kmer_positions[*end].0));
                chunks.push((*start, *end + 1));
            }
        }
        eprintln!("ok sangjin told me to contig {} kmer_positions.len() {}, chunk_positions {:?}, chunks {:?}", contig, kmer_positions.len(), chunk_positions, chunks);

        for (chunkdex, (left, right)) in chunks.iter().enumerate() {
            let left_pos = chunk_positions[chunkdex].0 - chunk_positions[chunkdex].0.min(30);
            let length = assembly.contig_sizes.get(contig).unwrap();
            let right_pos = chunk_positions[chunkdex].1 + (length - chunk_positions[chunkdex].1).min(30);
            let contig_start_pos = left_pos;

            let chunk_name = vec![
                contig_name.to_string(),
                (chunkdex + 1).to_string(),
                left_pos.to_string(),
                right_pos.to_string(),
            ]
            .join("_");
            eprintln!("left-right {}-{}", left, right);
            for ldex in *left..*right {
                //0..center.clusters[0].center.len() {
                // output is semi-vcf contig\tpos\t.\tREF\tALT\tqual\tfilter\tinfo\tformat\tsample

                let (pos, kmer) = kmer_positions[ldex];
                let reference;
                let alternate;
                let flip;
                reference = kmers.kmers.get(&kmer).unwrap().to_string();
                alternate = kmers.kmers.get(&Kmers::pair(kmer)).unwrap().to_string();
                match allele(kmer) {
                    Allele::Ref => {
                        flip = false;
                    }
                    Allele::Alt => {
                        //reference = kmers.kmers.get(&Kmers::pair(kmer)).unwrap().to_string();
                        //alternate = kmers.kmers.get(&kmer).unwrap().to_string();
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
                    (pos - contig_start_pos).to_string(),
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

fn get_phase_block_consistencies(phase_blocks: &HashMap<usize, Vec<usize>>, putative_phasing: &Vec<Option<bool>>, 
    kmer_positions: &Vec<(usize, i32)>, hic_mols: &Mols, hic_kmer_mols: &KmerMols) -> HashMap<(usize, usize), [u32; 4]> {
    let mut phase_block_consistencies: HashMap<(usize, usize), [u32; 4]> = HashMap::new();
    let mut blocks: Vec<usize> = Vec::new();
    
    
    let mut block_kmer_phasings: HashMap<usize, HashMap<i32, bool>> = HashMap::new();
    for (block_id, indices) in phase_blocks.iter() {
        blocks.push(*block_id);
        let kmer_phasings = block_kmer_phasings.entry(*block_id).or_insert(HashMap::new());
        for index in indices {
            let (_pos, kmer) = kmer_positions[*index];
            let canonical_kmer = Kmers::canonical_pair(kmer);
            if let Some(phase) = putative_phasing[*index] {
                kmer_phasings.insert(canonical_kmer, phase);
            }
        }
    }
    blocks.sort_by(|a, b| b.cmp(a));

   

    for phase1_blockdex in 0..blocks.len() {
        let phase_block1 = blocks[phase1_blockdex];
        //let (start1, end1) = phase_blocks.get(&phase_block1).unwrap();

        let block1_phasing = block_kmer_phasings.get(&phase_block1).unwrap();
        for phase2_blockdex in (phase1_blockdex + 1)..blocks.len() {
            let phase_block2 = blocks[phase2_blockdex];
            //let (start2, end2) = phase_blocks[phase_block2];
            let block2_phasing = block_kmer_phasings.get(&phase_block2).unwrap();
            for (kmer1, phase1) in block1_phasing.iter() {
                let phase1 = *phase1;
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

fn increment_consistency_counts(phase: bool, allele: i32, counts: &mut [u32;4]) {
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


fn add_kmer_and_update_phasing_consistency_counts(kmer_phasing_consistency_counts: &mut HashMap<i32, [u32;4]>, 
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

fn is_phasing_consistent(counts: &[u32;4], thresholds: &PhasingConsistencyThresholds, debug: bool) -> PhasingConsistency {
    let cis = (counts[0] + counts[1]) as f32; // 3
    let trans = (counts[2] + counts[3]) as f32; //270
    let total = cis + trans; // 273
    if debug {
        eprintln!("counts {:?}, cis {}, trans {}, total {}", counts, cis, trans, total);
    }
    if cis > trans { // false
        let consistent_percentage = cis/total;
        if debug {
            eprintln!("consistent percentage cis {} = {}/{}", consistent_percentage, trans, total);
        }
        if total > thresholds.min_count as f32 && consistent_percentage > thresholds.min_percent {
            let minor_allele = counts[0].min(counts[1]) as f32;
            if debug {
                eprintln!("minor_allele {} = {}.min({})", minor_allele, counts[2], counts[3]);
            }
            if minor_allele/cis > thresholds.minor_allele_fraction {
                if debug {
                    eprintln!("return is consistent true, cis true")
                }
                return PhasingConsistency{is_consistent: true, cis: true};
            }
        } 
    } else { // yes
        let consistent_percentage = trans/total; // 270/273 = 0.989
        if debug {
            eprintln!("consistent percentage trans {} = {}/{}", consistent_percentage, trans, total);
        }
        if total > thresholds.min_count as f32 && consistent_percentage > thresholds.min_percent { // 273 > 20 && 0.989 > 0.75
            let minor_allele = counts[2].min(counts[3]) as f32; // 120.min(150) == 120.0
            if debug {
                eprintln!("minor_allele {} = {}.min({})", minor_allele, counts[2], counts[3]);
            }
            if minor_allele/trans > thresholds.minor_allele_fraction { // 120.0 / 270.0 == .44
                if debug {
                    eprintln!("return is consistent true, cis false")
                }
                return PhasingConsistency{is_consistent: true, cis: false}; // return is_consistent true, cis false
            }
        } 
    }
    if debug {
        eprintln!("return is consistent false cis false");
    }
    
    PhasingConsistency{is_consistent: false, cis: false}
}



fn detect_sex_contigs(assembly: &Assembly, ccs_mols: &Mols, params: &Params, kmer_info: &Kmers) -> HashSet<i32> {
    let mut sex_contigs: HashSet<i32> = HashSet::new();
    let mut densities: Vec<(f32, f32, usize, usize, usize, usize, usize )> = Vec::new();
    let mut cov_sum: f32 = 0.0;
    let mut denom: f32 = 0.0;
    let mut density_sum: f32 = 0.0;
    eprintln!("starting pairwise");
    let pairwise_consistencies: HashMap<(i32, i32), [u32;4]> = get_pairwise_consistencies(&ccs_mols, assembly, true);
    eprintln!("pairwise done with {} entries", pairwise_consistencies.len());
    let thresholds = PhasingConsistencyThresholds {
        min_count: params.min_phasing_consistency_counts,
        min_percent: params.min_phasing_consistency_percent,
        minor_allele_fraction: params.min_minor_allele_fraction,
    };
    //eprintln!("ok how many contigs are there in the assembly {}", assembly.molecules.len());

    for contig_id in 1..(assembly.contig_names.len()) {
        let size = assembly.contig_sizes.get(&(contig_id as i32)).expect("I am actually going crazy");
        let kmers = match assembly.contig_kmers.get(&(contig_id as i32)) {
            Some(x) => x.len(),
            None => 0,
        };
        let mut consistent = 0;
        let mut inconsistent = 0;
        let mut consistent_kmers = 0;
        let mut inconsistent_kmers = 0;
        let kmer_positions = assembly.contig_kmers.get(&(contig_id as i32)).unwrap();
        for index1 in 0..kmer_positions.len() {
            let (pos1, kmer1) = kmer_positions[index1];
            let mut kmer1_consistent = 0.0;
            let mut kmer1_inconsistent = 0.0;
            let mut start = 0;
            if index1 > 100 { start = index1 - 100; }
            //eprintln!("ok index {} now we check {}-{}", index1, start, (index1 + 100).min(kmer_positions.len()));
            for index2 in start..(index1 + 100).min(kmer_positions.len()) {
                if index1 == index2 { continue; }
                let (pos2, kmer2) = kmer_positions[index2];
                //eprintln!("pos1 {} pos2 {}, {}-{} = {}", pos1, pos2,pos1.max(pos2), pos1.min(pos2), pos1.max(pos2)- pos1.min(pos2));
                if pos1.max(pos2) - pos1.min(pos2) < 5000 {
                    let kmer1 = Kmers::canonical_pair(kmer1.abs());
                    let kmer2 = Kmers::canonical_pair(kmer2.abs());
                    let key1 = kmer1.min(kmer2);
                    let key2 = kmer1.max(kmer2);
                    if let Some(count) = pairwise_consistencies.get(&(key1, key2)) {
                        let consistency = is_phasing_consistent(count, &thresholds, false);
                        let mut text = "NOT consistent";
                        if consistency.is_consistent { text = "IS consistent"; }
                        eprintln!("\t{}-{} = {:?} {}",pos1, pos2, count, text);
                        if consistency.is_consistent { kmer1_consistent += 1.0; consistent += 1; } else { kmer1_inconsistent += 1.0; inconsistent += 1; }
                    } else {
                        //eprintln!("no counts??? positions {}-{} length {}, kmers {} and {}", pos1, pos2, pos1.max(pos2) - pos1.min(pos2), 
                        //    kmer_info.kmers.get(&kmer1.abs()).unwrap(), kmer_info.kmers.get(&kmer2.abs()).unwrap());
                        kmer1_inconsistent += 1.0; inconsistent += 1;
                    }

                    
                } 
            }
            if kmer1_consistent + kmer1_inconsistent == 0.0 || kmer1_consistent/(kmer1_consistent + kmer1_inconsistent) > 0.25 { consistent_kmers += 1; } else { inconsistent_kmers += 1; }
            //eprintln!("contig {} kmer index {} position {} consistent with {} and inconsistent with {} so {}%", contig_id, index1, pos1, kmer1_consistent, kmer1_inconsistent, kmer1_consistent/(kmer1_consistent+kmer1_inconsistent));

        }

        //eprintln!("contig_id {}",contig_id);
        densities.push((params.contig_kmer_cov[contig_id], (kmers as f32)/(*size as f32), contig_id, consistent, inconsistent, consistent_kmers, inconsistent_kmers));
        let size = *size as f32;
        cov_sum += params.contig_kmer_cov[contig_id] * size;
        denom += size;
        density_sum += kmers as f32;
    }

    let avg_cov = cov_sum / denom;
    let avg_density = density_sum / denom;

    eprintln!("detecting sex contigs. mean kmer count is {} and mean paired kmer density is {}", avg_cov, avg_density);
    eprintln!("kmer_depth\thet_kmer_density\tcontig_id\tcontig_name\tcontig_length\tcontig_classification\tconsistent_links\tinconsistent_links\tconsistent_kmers\tinconsistent_kmers\tsex_contig_cov_cutoff\tsex_density_cutoff");
    for (depth, density, contig, consistent, inconsistent, consistent_kmers, inconsistent_kmers) in densities.iter() {
        let length = *assembly.contig_sizes.get(&(*contig as i32)).unwrap();
        let mut sex = "autosome";
        if length > params.min_contig_length && *depth < params.sex_contig_cov_cutoff * avg_cov 
            && *density < params.sex_contig_het_kmer_density_cutoff * avg_density {
                sex_contigs.insert(*contig as i32);
            sex = "sex";
        }
        eprintln!("{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}", depth, density, contig, 
            assembly.contig_names[*contig as usize], 
            length, sex, consistent, inconsistent, consistent_kmers, inconsistent_kmers, params.sex_contig_cov_cutoff * avg_cov, 
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

    let min_minor_allele_fraction = params.value_of("min_minor_allele_fraction").unwrap_or("0.25");
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