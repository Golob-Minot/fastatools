#!/usr/bin/env python
import argparse
import fastalite
import logging

# As the name implies, given a set of fasta / fastq files (min 2), combine them into one fasta file.
#  Check at least to be sure no overlapping IDs. Optionally check sequences themselves.
#  Optionally handle paired reads (being sure to include each pair)


def main():
    args_parser = argparse.ArgumentParser(
        description="""Given a set of fasta / fastq files (min 2), combine them into one fasta file.
        Check at least to be sure no overlapping IDs. Optionally check for overlapping sequences too.
        """
    )

    args_parser.add_argument(
        'files',
        help='Files to be combined. Minimum of 2 required',
        nargs='+',
        type=fastalite.Opener(mode='r')
    )
    args_parser.add_argument(
        '--fastq',
        '-q',
        help='input and outputs are fastq (not fasta)',
        action='store_true'
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
        help='File into which we should place our combined reads',
        required=True,
        type=fastalite.Opener(mode='w')
    )

    args = args_parser.parse_args()
    logging.basicConfig(level=logging.INFO)

    if len(args.files) < 2:
        logging.error("Only one file given. Nothing to do.")
        return -1

    out_h = args.output

    if args.check_seq:
        seq_ids = set()
        seqs = set()
        for file_h in args.files:
            if args.fastq:
                reader = fastalite.fastqlite(file_h)
            else:
                reader = fastalite.fastalite(file_h)
            for sr in reader:
                if sr.id not in seq_ids and sr.seq not in seqs:
                    seq_ids.add(sr.id)
                    seqs.add(sr.seq)
                    if args.fastq:
                        out_h.write("@{} {}\n{}\n+\n{}\n".format(
                            sr.id,
                            sr.description,
                            sr.seq,
                            sr.qual
                        ))
                    else:
                        out_h.write(">%s %s\n%s\n" % (sr.id, sr.description, sr.seq))

    else:  # just IDs
        seq_ids = set()
        for file_h in args.files:
            if args.fastq:
                reader = fastalite.fastqlite(file_h)
            else:
                reader = fastalite.fastalite(file_h)
            for sr in reader:
                if sr.id not in seq_ids:
                    seq_ids.add(sr.id)
                    if args.fastq:
                        out_h.write("@{} {}\n{}\n+\n{}\n".format(
                            sr.id,
                            sr.description,
                            sr.seq,
                            sr.qual
                        ))
                    else:
                        out_h.write(">%s %s\n%s\n" % (sr.id, sr.description, sr.seq))

if __name__ == "__main__":
    main()
