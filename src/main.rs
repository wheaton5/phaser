#[macro_use]
extern crate clap;
extern crate bio;
extern crate disjoint_set;
extern crate hashbrown;
extern crate phasst_lib;

use bio::io::fasta;
use bio::io::fasta::Record;
use bio::utils::TextSlice;
use std::path::Path;

use phasst_lib::{
    load_assembly_kmers, load_hic, load_hifi, load_linked_read_barcodes, Assembly, HicMols,
    HifiMols, Kmers, LinkedReadBarcodes,
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

}

fn detect_sex_contigs(assembly: &Assembly, params: &Params) -> HashSet<i32> {
    let mut sex_contigs: HashSet<i32> = HashSet::new();
    let mut densities: Vec<(f32, f32, i32)> = Vec::new();
    let mut cov_sum: f32 = 0.0;
    let mut denom: f32 = 0.0;
    let mut density_sum: f32 = 0.0;


    for (contig, kmers) in assembly.molecules.iter() {
        let size = assembly.contig_sizes.get(contig).expect("I am actually going crazy");
        densities.push((params.contig_kmer_cov[*contig as usize], (kmers.len() as f32)/(*size as f32), *contig));
        let size = *size as f32;
        cov_sum += params.contig_kmer_cov[*contig as usize] * size;
        denom += size;
        density_sum += kmers.len() as f32;
    }

    let avg_cov = cov_sum / denom;
    let avg_density = density_sum / denom;

    //densities.sort_by(|a, b| a.partial_cmp(b).unwrap());
    for (depth, density, contig) in densities.iter() {
        let mut sex = "autosomal";
        if *depth < params.sex_contig_cov_cutoff * avg_cov 
            && *density < params.sex_contig_het_kmer_density_cutoff * avg_density {
                sex_contigs.insert(*contig);
            sex = "sex";
        }
        eprintln!("{}\t{}\t{}\t{}\t{}\t{}", depth, density, contig, assembly.contig_names[*contig as usize], assembly.contig_sizes.get(contig).unwrap(), sex);

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
    break_window: usize,
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
        ploidy: ploidy,
        min_hic_links: min_hic_links,
        break_window: break_window,
    }
}