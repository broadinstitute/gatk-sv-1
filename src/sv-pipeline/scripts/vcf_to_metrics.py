#!/bin/python

from pysam import VariantFile

import argparse
import sys
from typing import Optional, List, Text, TextIO

# Constants
NA_VALUE = "NaN"
PE_HEADER = "\t".join(["name", "PEQ", "PECS"])
SR_HEADER = "\t".join(["name", "coord", "pos", "SRQ", "SRCS"])
PESR_HEADER = "\t".join(["name", "PESRQ", "PESRCS"])
BAF_HEADER = "\t".join(["name", "BAFDEL", "BAFDUP"])


# Helper functions to write metric output rows for each type of evidence
def _write_sr_line(f: TextIO, record: VariantFile) -> None:
  if 'SRQ' not in record.info:
    return
  if 'SR1POS' in record.info and record.info['SR1POS'] is not None:
    sr1pos = record.info['SR1POS']
  else:
    sr1pos = record.pos
  if 'SR2POS' in record.info and record.info['SR2POS'] is not None:
    sr2pos = record.info['SR2POS']
  else:
    sr2pos = record.stop
  cols = [record.id, "posA", sr1pos, record.info['SR1Q'], record.info['SR1CS']]
  f.write("\t".join(str(c) for c in cols) + "\n")
  cols = [record.id, "posB", sr2pos, record.info['SR2Q'], record.info['SR2CS']]
  f.write("\t".join(str(c) for c in cols) + "\n")
  cols = [record.id, "sum", 0, record.info['SRQ'], record.info['SRCS']]
  f.write("\t".join(str(c) for c in cols) + "\n")


def _write_pe_line(f: TextIO, record: VariantFile) -> None:
  if 'PEQ' not in record.info:
    return
  cols = [record.id, record.info['PEQ'], record.info['PECS']]
  f.write("\t".join(str(c) for c in cols) + "\n")


def _write_pesr_line(f: TextIO, record: VariantFile) -> None:
  if 'PESRQ' not in record.info:
    return
  cols = [record.id, record.info['PESRCS'], record.info['PESRQ']]
  f.write("\t".join(str(c) for c in cols) + "\n")


def _write_baf_line(f: TextIO, record: VariantFile) -> None:
  if 'BAFDEL' in record.info:
    cols = [record.id, record.info['BAFDEL'], NA_VALUE]
    f.write("\t".join(str(c) for c in cols) + "\n")
  elif 'BAFDUP' in record.info:
    cols = [record.id, NA_VALUE, record.info['BAFDUP']]
    f.write("\t".join(str(c) for c in cols) + "\n")


# Argument definitions
def __parse_arguments(argv: List[Text]) -> argparse.Namespace:
  parser = argparse.ArgumentParser(
    description="Extracts SV metric annotations from a VCF and writes them to tsv files",
    formatter_class=argparse.ArgumentDefaultsHelpFormatter
  )
  parser.add_argument("--vcf", type=str, required=True,
                      help="Input vcf path")
  parser.add_argument("--pe-out", type=str, required=True,
                      help="PE metrics output path")
  parser.add_argument("--sr-out", type=str, required=True,
                      help="SR metrics output path")
  parser.add_argument("--pesr-out", type=str, required=True,
                      help="PESR metrics output path")
  parser.add_argument("--baf-out", type=str, required=True,
                      help="BAF metrics output path")
  if len(argv) <= 1:
    parser.parse_args(["--help"])
    sys.exit(0)
  parsed_arguments = parser.parse_args(argv[1:])
  return parsed_arguments


def main(argv: Optional[List[Text]] = None):
  if argv is None:
    argv = sys.argv
  args = __parse_arguments(argv)
  with VariantFile(args.vcf) as vcf, \
          open(args.pe_out, 'w') as fpe, \
          open(args.sr_out, 'w') as fsr, \
          open(args.pesr_out, 'w') as fpesr, \
          open(args.baf_out, 'w') as fbaf:
    fpe.write(PE_HEADER + "\n")
    fsr.write(SR_HEADER + "\n")
    fpesr.write(PESR_HEADER + "\n")
    fbaf.write(BAF_HEADER + "\n")
    for r in vcf:
      _write_pe_line(fpe, r)
      _write_sr_line(fsr, r)
      _write_pesr_line(fpesr, r)
      _write_baf_line(fbaf, r)


if __name__ == "__main__":
  main()
