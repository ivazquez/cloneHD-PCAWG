#!/usr/bin/env python2
from __future__ import print_function

# Requires PyVCF. To install: pip2 install pyvcf
import vcf
import argparse
import csv
from collections import defaultdict
import random
import sys

class VariantParser(object):
  def __init__(self):
    self._vcf_filename = None

  def list_variants(self):
    variants = self._filter(self._vcf_filename)
    variants_and_reads = []
    for variant in variants:
      ref_reads, alt_reads, total_reads = self._calc_read_counts(variant)
      variants_and_reads.append((variant, ref_reads, alt_reads, total_reads))
    return variants_and_reads

  def _parse_vcf(self, vcf_filename):
    vcfr = vcf.Reader(filename=vcf_filename)
    records = []
    for variant in vcfr:
      variant.CHROM = variant.CHROM.lower()
      # some VCF dialects prepend "chr", some don't
      # remove the prefix to standardize
      if variant.CHROM.startswith('chr'):
        variant.CHROM = variant.CHROM[3:]
      records.append(variant)
    return records

  def _is_good_chrom(self, chrom):
    if chrom in [str(i) for i in range(1, 23)]:
      return True
    elif chrom in ['x', 'y']:
      return True
    else:
      return False

  def _does_variant_pass_filters(self, variant):
    if variant.FILTER is None:
      return True
    if len(variant.FILTER) > 0:
      # variant failed one or more filters
      return False
    return True

  def _filter(self, vcf_filename):
    variants = []

    all_variants = self._parse_vcf(vcf_filename)

    for variant in all_variants:
      if not self._is_good_chrom(variant.CHROM):
        continue
      if not self._does_variant_pass_filters(variant):
        continue
      variants.append(variant)
    return variants

  def _get_sample_index(self, variant, sample=None):
    """Find the index of the sample.
    """
    if self._sample:
      sample_is = [i for i, s in enumerate(variant.samples) if s.sample == sample]
      assert len(sample_is) == 1, "Did not find sample name %s in samples" % sample
      return sample_is[0]
    else:
      # don't make this -1, as some code assumes it will be >= 0
      return len(variant.samples) - 1

class MutectSmchetParser(VariantParser):
  def __init__(self, vcf_filename, sample=None):
    self._vcf_filename = vcf_filename
    self._sample = sample

  def _calc_read_counts(self, variant):
    sample_i = self._get_sample_index(variant, self._sample)
    ref_reads = int(variant.samples[sample_i]['AD'][0])
    alt_reads = int(variant.samples[sample_i]['AD'][1])
    total_reads = ref_reads + alt_reads

    return (ref_reads, alt_reads, total_reads)

class VariantFormatter(object):
  def __init__(self):
    self._counter = 0

  def _split_types(self, genotype):
    types = [int(e) for e in genotype.split('/')]
    if len(types) != 2:
      raise Exception('Not diploid: %s' % types)
    return types

  def format_variants(self, variant_list):
    variant_list.sort(key = lambda v: variant_key(v[0]))
    
    for variant, ref_reads, alt_reads, total_reads in variant_list:
      snv_id = 's%s' % self._counter
      chrom, pos = variant_key(variant)

      yield {
        'chrom': chrom,
        'pos': pos,
        'ref_reads': ref_reads,
        'alt_reads': alt_reads,
        'total_reads': total_reads
      }
      self._counter += 1
      
  def write_variants(self, variants, outfn):
    with open(outfn, 'w') as outf:
      for variant in variants:
        vals = (
          'chrom',
          'pos',
          'alt_reads',
          'total_reads',
        )
        vals = [variant[k] for k in vals]
        print('\t'.join([str(v) for v in vals]), file=outf)

def restricted_float(x):
  x = float(x)
  if x < 0.0 or x > 1.0:
    raise argparse.ArgumentTypeError('%r not in range [0.0, 1.0]' % x)
  return x

def variant_key(var):
  chrom = var.CHROM
  if chrom == 'x':
    chrom = 23
  elif chrom == 'y':
    chrom = 24
  elif chrom.isdigit():
    chrom = int(chrom)
  else:
    chrom = 999
  return (chrom, var.POS)

def main():
  parser = argparse.ArgumentParser(
    description='Create snv_data.txt input files for cloneHD from VCF data.',
    formatter_class=argparse.ArgumentDefaultsHelpFormatter
  )
  parser.add_argument('-i', '--vcf', dest='vcf_file', required=True, 
    help='Input VCF file')
  parser.add_argument('-t', '--vcf-type', dest='input_type', required=True, choices=('mutect-smchet',),
    help='Type of VCF file')
  parser.add_argument('-s', '--sample', dest='sample',
    help='Name of the sample in the input VCF file. Defaults to last sample if not specified.')
  parser.add_argument('-o', '--snv', dest='snv_file', default='snv_data.txt',
    help='Output SNV file')
  args = parser.parse_args()
  
  if args.input_type == 'mutect-smchet':
    parser = MutectSmchetParser(args.vcf_file, args.sample)
  variant_list = parser.list_variants()
  
  formatter = VariantFormatter()
  variant_format = list(formatter.format_variants(variant_list))
  
  formatter.write_variants(variant_format, args.snv_file)

if __name__ == '__main__':
  main()