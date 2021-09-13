#!/usr/bin/env python

import argparse
from collections import Counter
import gzip
import pysam
import sys


def main():
    parser = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument('revise_vcf_lines', type=argparse.FileType('r'))
    parser.add_argument('normal_revise_vcf')
    parser.add_argument('famfile', type=argparse.FileType('r'))
    parser.add_argument('sexchr_revise')
    parser.add_argument('multi_geno_ids_txt')
    parser.add_argument('outlier_samples_list', type=argparse.FileType('r'))
    args = parser.parse_args()

    # load the revised lines and index by ID
    revised_lines_by_id = {}
    with pysam.VariantFile(args.revise_vcf_lines) as revise_vcf:
        header2 = revise_vcf.header
        for record in revise_vcf:
            revised_lines_by_id[record.id] = record
    print("loaded {} revised lines".format(len(revised_lines_by_id)), file=sys.stderr)

    outlier_samples = set([line.rstrip() for line in args.outlier_samples_list if not line.isspace()])
    print("loaded {} outlier samples".format(len(outlier_samples)), file=sys.stderr)

    male_samples = set()
    for line in args.famfile:
        if line.isspace():
            continue
        fields = line.rstrip().split("\t")
        if fields[4] == '1':
            male_samples.add(fields[1])
    print("identified {} male samples".format(len(male_samples)), file=sys.stderr)

    if args.sexchr_revise.endswith(".gz"):
        sexchr_revise = set([line.rstrip() for line in gzip.open(args.sexchr_revise, 'rt')])
    else:
        sexchr_revise = set([line.rstrip() for line in open(args.sexchr_revise, 'rt')])
    print("{} sites to revise on sex chromosomes".format(len(sexchr_revise)), file=sys.stderr)

    if args.multi_geno_ids_txt.endswith(".gz"):
        multi_geno_ids = set([line.rstrip() for line in gzip.open(args.multi_geno_ids_txt, 'rt')])
    else:
        multi_geno_ids = set([line.rstrip() for line in open(args.multi_geno_ids_txt, 'rt')])
    print("{} multiallelic sites".format(len(multi_geno_ids)), file=sys.stderr)

    normal_vcf = pysam.VariantFile(args.normal_revise_vcf)

    NEW_HEADER_LINES = ['##ALT=<ID=CNV,Description="Copy Number Polymorphism">',
                        '##FORMAT=<ID=CNQ,Number=1,Type=Integer,Description="Read-depth genotype quality">',
                        '##FORMAT=<ID=CN,Number=1,Type=Integer,Description="Predicted copy state">',
                        '##FILTER=<ID=PESR_GT_OVERDISPERSION,Description="High PESR dispersion count">',
                        '##FILTER=<ID=MULTIALLELIC,Description="Multiallelic site">']

    # # Add metadata lines for annotations
    header1 = normal_vcf.header

    for f in NEW_HEADER_LINES:
        header1.add_line(f)
        header2.add_line(f)

    non_outlier_samples = set([s for s in header1.samples if not s in outlier_samples])
    vf_1 = max(len(non_outlier_samples) * 0.01, 2)

    biallelic_gts = { (1, 1), (0, 0), (0, 1), (None, None)}

    print("reformatting records", file=sys.stderr)
    cleangq_filename = "cleanGQ.vcf.gz"
    cleanqg_out = pysam.VariantFile(cleangq_filename, 'w', header=normal_vcf.header)

    mulitallelic_filename = "multiallelic.vcf.gz"
    multiallelic_out = pysam.VariantFile(mulitallelic_filename, 'w', header=normal_vcf.header)

    for idx, record in enumerate(normal_vcf):
        multi_del = False
        multi_dup = False
        gt4_copystate = False
        gt5kb_dup = False
        gt5kb_del = False
        if (idx - 1) % 1000 == 0:
            print("processed {} records".format(idx), file=sys.stderr)
        if record.id in revised_lines_by_id:
            record = revised_lines_by_id[record.id]
        if record.info.get('SVTYPE', None) == 'DEL':
            if abs(record.stop - record.pos) >= 1000:
                sample_cn_map = {s: record.samples[s]['RD_CN'] for s in non_outlier_samples}
                if len([s for s in sample_cn_map if (sample_cn_map[s] is not None and sample_cn_map[s] > 3)]) > vf_1:
                    multi_del = True
            gts = [record.samples[s]['GT'] for s in non_outlier_samples]
            if len([gt for gt in gts if gt not in biallelic_gts]) > 0:
                gt5kb_del = True
            if abs(record.stop - record.pos) >= 5000:
                if not multi_del:
                    gt5kb_del = True

        if record.info.get('SVTYPE', None) == 'DUP':
            if abs(record.stop - record.pos) >= 1000:
                sample_cn_map = {s: record.samples[s]['RD_CN'] for s in non_outlier_samples}
                if len([s for s in sample_cn_map if (sample_cn_map[s] is not None and sample_cn_map[s] > 4)]) > vf_1:
                    multi_dup = True
                if len([x for x in Counter(sample_cn_map.values()) if (x is not None and (x < 1 or x > 4))]) > 4:
                    gt4_copystate
                if len([s for s in sample_cn_map if (sample_cn_map[s] is not None and (sample_cn_map[s] < 1 or sample_cn_map[s] > 4) and
                                                     gt4_copystate)]) > vf_1:
                    multi_dup = True
            gts = [record.samples[s]['GT'] for s in non_outlier_samples]
            if len([gt for gt in gts if gt not in biallelic_gts]) > 0:
                gt5kb_dup = True
            if abs(record.stop - record.pos) >= 5000:
                if not multi_dup:
                    gt5kb_dup = True

        if gt5kb_del:
            for sample in record.samples:
                if not record.samples[sample]['GQ'] is None and record.samples[sample]['RD_CN'] >= 2:
                    record.samples[sample]['GT'] = (0, 0)
                elif not record.samples[sample]['GQ'] is None and record.samples[sample]['RD_CN'] == 1:
                    record.samples[sample]['GT'] = (0, 1)
                elif not record.samples[sample]['GQ'] is None:
                    record.samples[sample]['GT'] = (1, 1)  # RD_CN 0 DEL

        if gt5kb_dup:
            for sample in record.samples:
                if not record.samples[sample]['GQ'] is None and record.samples[sample]['RD_CN'] <= 2:
                    record.samples[sample]['GT'] = (0, 0)
                elif not record.samples[sample]['GQ'] is None and record.samples[sample]['RD_CN'] == 3:
                    record.samples[sample]['GT'] = (0, 1)
                elif not record.samples[sample]['GQ'] is None:
                    record.samples[sample]['GT'] = (1, 1)  # RD_CN > 3 DUP

        if record.id in multi_geno_ids:
            record.filter.add('PESR_GT_OVERDISPERSION')

        if multi_del or multi_dup:
            record.filter.add('MULTIALLELIC')
            for j, sample in enumerate(record.samples):
                record.samples[sample]['GT'] = None
                record.samples[sample]['GQ'] = None
                record.samples[sample]['CN'] = record.samples[sample]['RD_CN']
                record.samples[sample]['CNQ'] = record.samples[sample]['RD_GQ']

        if len(record.filter) > 1 and 'PASS' in record.filter:
            del record.filter['PASS']

        if 'MULTIALLELIC' in record.filter and ('<DUP>' in record.alts or '<DEL>' in record.alts):
            record.alts = ('<CNV>',)
            record.info['SVTYPE'] = 'CNV'

        if record.id in sexchr_revise:
            for sample in record.samples:
                if sample in male_samples:
                    cn = int(record.samples[sample]['RD_CN'])
                    if cn is not None and cn > 0:
                        record.samples[sample]['RD_CN'] = cn - 1
                        if 'CN' in record.samples[sample]:
                           record.samples[sample]['CN'] = cn - 1  # the old script didn't do this but I think it should

        cleanqg_out.write(record)

        if 'MULTIALLELIC' in record.filter:
            multiallelic_out.write(record)

    cleanqg_out.close()
    multiallelic_out.close()
    print("done", file=sys.stderr)


if __name__ == '__main__':
    main()