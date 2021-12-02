#!/usr/bin/env python

import argparse
import csv
import lzma
import time

import pandas as pd
from path import Path

from Bio import AlignIO, SeqIO

# import bjorn_support as bs
import mutations as bm
import data as bd


''' Hacked up from the version in bjorn_support to do xz.
'''
def open_fasta_xz(fasta_filepath, is_aligned=False):
  handle = lzma.open(fasta_filepath, "rt")
  if is_aligned:
    cns = AlignIO.read(handle, 'fasta')
  else:
    cns = SeqIO.parse(handle, 'fasta')
  # for record in cns:
  #   print(record.id)
  return handle, cns


''' Modified to include insertions, and to write TSV instead of CSV. [fry 210713]
'''
def process_mutations(alignment_filepath, patient_zero, output_path, data_src='gisaid_feed'):
  print(f"Loading alignment file at {alignment_filepath}")
  t0 = time.time()
  # msa_data = bs.load_fasta(alignment_filepath, is_aligned=True, is_gzip=False)
  msa_stream, msa_data = open_fasta_xz(alignment_filepath, is_aligned=True)
  msa_load_time = time.time() - t0

  print("Identifying insertions...")
  t0 = time.time()
  inserts, _ = bm.identify_insertions_per_sample(msa_data,
                                                 #gisaid_meta,
                                                 gene2pos=bd.GENE2POS,
                                                 data_src=data_src,
                                                 patient_zero=patient_zero
                                                 #test=is_test
                                                 )
  inserts_time = time.time() - t0

  print("Identifying substitution-based mutations...")
  t0 = time.time()
  subs, _ = bm.identify_replacements_per_sample(msa_data,
                                                #gisaid_meta,
                                                gene2pos=bd.GENE2POS,
                                                data_src=data_src,
                                                min_seq_len=20000,
                                                patient_zero=patient_zero
                                                #test=is_test
                                                )
  subs_time = time.time() - t0

  print("Identifying deletion-based mutations...")
  t0 = time.time()
  dels, _ = bm.identify_deletions_per_sample(msa_data,
                                             #gisaid_meta,
                                             gene2pos=bd.GENE2POS,
                                             data_src=data_src,
                                             min_del_len=1,
                                             max_del_len=500,
                                             min_seq_len=20000,
                                             patient_zero=patient_zero
                                             #test=is_test
                                             )
  dels_time = time.time() - t0

  # QC FILTER: remove seqs with >500 nt deletions
  # dels = dels.loc[dels['del_positions'].str.len()<500]

  print(inserts.shape)
  print(subs.shape)
  print(dels.shape)

  # muts = pd.concat([subs, dels])
  muts = pd.concat([inserts, subs, dels])

  # create the is_synonymous column by comparing ref_aa with alt_aa
  muts['is_synonymous'] = False
  muts.loc[muts['ref_aa']==muts['alt_aa'], 'is_synonymous'] = True

  # make sure we have the expected columns, because if there are
  # no insertions (or deletions, or substitutions for that matter),
  # those columns will be missing
  expected_columns = [
    'absolute_coords', 'alt_aa', 'alt_codon', 'codon_num', 'codon_start',
    'del_len', 'del_positions', 'del_seq', 'deletion_codon_coords',
    'gene', 'gene_start_pos', 'idx', 'ins_len', 'ins_positions',
    'is_frameshift', 'is_synonymous', 'mutation', 'next_5nts',
    'pos', 'pos_in_codon', 'prev_5nts', 'ref_aa', 'ref_codon',
    'relative_coords', 'replacements', 'seq_len', 'type'
  ]
  missing_columns = set(expected_columns)
  for col in muts.columns:
    missing_columns.remove(col)

  # add empty columns for the missing fellers
  for missing in missing_columns:
    print(f'adding empty column for {missing}')
    muts[missing] = ''

  print(muts.shape)
  # muts = muts.astype(str) TAKES FOREVER
  # muts_filename = alignment_filepath.replace('.aligned.fasta', f'_{date}.mutations.csv')

  # put columns into alphabetical order, because the ordering seems to be arbitrary,
  # and this at least keeps them in some sort of comparable order between runs
  muts.sort_index(axis=1, inplace=True)
  # muts.to_csv(output_path, index=False)  # old version was CSV, so we'll keep using that; oh well
  with lzma.open(output_path, 'wt') as ofp:
    muts.to_csv(ofp, sep='\t', quoting=csv.QUOTE_NONE, index=False)  # TSV version

  print(f"Mutations extracted from {alignment_filepath} and saved in {output_path}\n")

  # TODO should be earlier; how soon can this be closed?
  msa_stream.close()


if __name__=="__main__":
  # COLLECTING USER PARAMETERS
  parser = argparse.ArgumentParser()
  parser.add_argument("-i", "--input", type=str, required=True, help="FASTA filepath containing aligned sequences")
  # parser.add_argument("-m", "--meta", type=str, required=True, help="Gzipped TSV filepath containing sequence metadata")
  parser.add_argument("-r", "--patient-zero", type=str, default="NC_045512.2", help="Sample name of reference sequence")
  parser.add_argument("-d", "--data-src", type=str, default="gisaid_feed", help="Data source")
  parser.add_argument("-o", "--outfp", type=str, required=True, help="Output filepath storing mutation information")
  args = parser.parse_args()
  alignment_filepath = args.input
  # gisaid_meta = args.meta
  patient_zero = args.patient_zero
  data_src = args.data_src
  out_fp = Path(args.outfp)

  process_mutations(alignment_filepath, patient_zero, data_src, out_fp)
