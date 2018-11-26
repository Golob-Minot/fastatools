#!/usr/bin/env python
import argparse
import fastalite
import logging

# Given two fasta files, return only those reads in A that are NOT in B.
# Minimally considers sequence IDs. Can optionally also consider the actual sequences


def main():
    args_parser = argparse.ArgumentParser(
        description="""Given two fasta files, return only those reads in A that are NOT in B.
        Minimally considers sequence IDs. Can optionally also consider the actual sequences.
        """
    )

    args_parser.add_argument(
        'fasta_A',
        help='Fasta file A',
        type=argparse.FileType(mode='r')
    )
    args_parser.add_argument(
        'fasta_B',
        help='Fasta file B',
        type=argparse.FileType(mode='r')
    )
    args_parser.add_argument(
        '--check-seq',
        '-s',
        help='Also check to be sure sequences are not repeated. Default is to only check for repeated IDs',
        action='store_true'
    )
    args_parser.add_argument(
        '--output',
        '-o',
        help='Output file (fasta)',
        required=True,
        type=argparse.FileType('w')
    )

    args = args_parser.parse_args()
    logging.basicConfig(level=logging.INFO)

    out_h = args.output

    if args.check_seq:
        seq_ids = set()
        seqs = set()
        for sr in fastalite.fastalite(args.fasta_B):
            seq_ids.add(sr.id)
            seqs.add(sr.seq)
        for sr in fastalite.fastalite(args.fasta_A):
            if sr.id not in seq_ids and sr.seq not in seqs:
                out_h.write(">%s %s\n%s\n" % (sr.id, sr.description, sr.seq))
    else:  # just IDs
        seq_ids = set()
        for sr in fastalite.fastalite(args.fasta_B):
            seq_ids.add(sr.id)
        for sr in fastalite.fastalite(args.fasta_A):
            if sr.id not in seq_ids:
                out_h.write(">%s %s\n%s\n" % (sr.id, sr.description, sr.seq))                    

if __name__ == "__main__":
    main()
