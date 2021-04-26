import pandas as pd
from path import Path
from shutil import copy, move
from Bio import SeqIO, AlignIO, Phylo
import argparse
import glob
import subprocess
from multiprocessing import Pool
from itertools import repeat
import os
from datetime import datetime as dt
import bjorn_support as bs
import mutations as bm
# from bjorn_support import concat_fasta, align_fasta, compute_tree, map_gene_to_pos, load_fasta
# from mutations import identify_replacements, identify_deletions, identify_insertions
import data as bd
import json


## MAIN

if __name__=="__main__":
    with open('config.json', 'r') as f:
        config = json.load(f)
    # out_dir = Path(config['release_outdir'])
    # date = config['biolabs_date']
    # meta_fp = config['biolabs_meta']
    parser = argparse.ArgumentParser()
    parser.add_argument("-o", "--out-dir",
                        type=str,
                        help="Path to folder where results are to be saved")
    parser.add_argument("-d", "--date",
                        type=str,
                        help="Date assigned to the sequencing run")
    parser.add_argument("-m", "--metadata",
                        type=str,
                        help="Path to file containing biolabs metadata")
    args = parser.parse_args()
    out_dir = Path(args.out_dir)
    date = args.date
    meta_fp = args.metadata
    num_cpus = config['num_cpus']
    ref_path = config['reference_filepath']
    patient_zero = config['patient_zero']
    min_coverage = config['min_coverage']
    # msa_fp = config['msa_filepath']
    nonconcerning_genes = config['nonconcerning_genes']
    nonconcerning_mutations = config['nonconcerning_mutations']
    if not Path.isdir(out_dir):
        Path.mkdir(out_dir);
    msa_dir = out_dir/'msa'
    if not Path.isdir(msa_dir):
        Path.mkdir(msa_dir);
    seqs_dir = Path(out_dir/'fa')
    accepted_seqs_dir = Path(out_dir/'fa_accepted')
    if not Path.isdir(accepted_seqs_dir):
        Path.mkdir(accepted_seqs_dir);
    intro = pd.read_excel(meta_fp, sheet_name='Instructions')
    cov = pd.read_excel(meta_fp, sheet_name='Coverage')
    meta = pd.read_excel(meta_fp, sheet_name='Submissions', skiprows=1)
    accepted_samples = cov.loc[cov['Uniformity of base cov.']>min_coverage, 'Biolab Trans. #'].tolist()
    meta['sample_id'] = meta['FASTA filename'].apply(lambda x : x.split('_')[1].split('.')[0]).astype(int)
    accepted_sample_filenames = meta.loc[meta['sample_id'].isin(accepted_samples), 'FASTA filename'].tolist()
    meta = meta.loc[meta['sample_id'].isin(accepted_samples)].drop(columns=['sample_id'])
    meta[['Gender', 'Patient age', 'Patient status']] = 'N/A'
    meta.to_csv(f'{out_dir}/raw_gisaid_metadata.csv', index=False)
    for sample_filename in accepted_sample_filenames:
        copy(f'{seqs_dir}/{sample_filename}', accepted_seqs_dir)
    copy(ref_path, accepted_seqs_dir)
    fasta_fps = glob.glob(f"{accepted_seqs_dir}/*.fasta")
    out_fasta_fp = f"{msa_dir}/{date}_release.fa"
    all_sequences = []
    ref_seq = SeqIO.read(ref_path, 'fasta')
    all_sequences.append(ref_seq)
    for fp in fasta_fps:
        rec = SeqIO.read(fp, 'fasta')
        all_sequences.append(rec)
    SeqIO.write(all_sequences, out_fasta_fp, 'fasta')
    # seqs_fp = bs.concat_fasta(accepted_seqs_dir, msa_dir/out_dir.basename());
    # load concatenated sequences
    # cns_seqs = SeqIO.parse(msa_dir/out_dir.basename()+'.fa', 'fasta')
    # cns_seqs = list(cns_seqs)
    # print(len(cns_seqs))
    msa_fp = out_fasta_fp.split('.')[0] + '_aligned.fa'
    if not Path.isfile(Path(msa_fp)):
        msa_fp = bs.align_fasta(out_fasta_fp, msa_fp, num_cpus=num_cpus);
    # load multiple sequence alignment
    msa_data = bs.load_fasta(msa_fp, is_aligned=True)
    # identify insertions
    insertions = bm.identify_insertions(msa_data, 
                                        # meta_fp=meta_fp, 
                                        patient_zero=patient_zero, 
                                        min_ins_len=1)
    # save insertion results to file
    insertions.to_csv(out_dir/'insertions.csv', index=False)
    # identify substitution mutations
    substitutions = bm.identify_replacements(msa_data,
                                # meta_fp=meta_fp,
                                patient_zero=patient_zero)
    # save substitution results to file
    substitutions.to_csv(out_dir/'replacements.csv', index=False)
    # identify deletions
    deletions = bm.identify_deletions(msa_data,
                                    # meta_fp=meta_fp,
                                    patient_zero=patient_zero,
                                    min_del_len=1)
    # save deletion results to file
    deletions.to_csv(out_dir/'deletions.csv', index=False)
    # identify samples with suspicious INDELs and/or substitutions
    sus_ids, sus_muts = bm.identify_samples_with_suspicious_mutations(substitutions, 
                                                                        deletions, 
                                                                        pd.DataFrame(),
                                                                        nonconcerning_genes,
                                                                        nonconcerning_mutations)
    sus_muts.to_csv(out_dir/'suspicious_mutations.csv', index=False)
    print(msa_fp)
