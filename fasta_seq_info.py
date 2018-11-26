#!/usr/bin/env python
import argparse
import fastalite
import csv
import logging
import sys

# Given a FASTA file(s) and sequence information csv file(s),
# filter the sequence information to only include rows for reads in the fasta file(s).


def main():
    args_parser = argparse.ArgumentParser(
        description="""Given a FASTA file(s) and sequence information csv file(s)
        filter the sequence information to only include rows for reads in the fasta file(s).
        """
    )

    args_parser.add_argument(
        'fasta',
        nargs='+',
        help='Fasta file(s)',
        type=argparse.FileType(mode='r')
    )
    args_parser.add_argument(
        '--sequence-info',
        '-si',
        nargs='+',
        help='Sequence information file(s) (csv format)',
        type=argparse.FileType(mode='r'),
        required=True
    )
    args_parser.add_argument(
        '--output',
        '-o',
        help='File into which we should place our filtered sequence information',
        required=True,
        type=argparse.FileType('w')
    )

    args = args_parser.parse_args()
    logging.basicConfig(level=logging.INFO)

    #  Create readers for each of our incoming seq_info files
    seq_info_readers = [
        csv.DictReader(seq_info_h) for seq_info_h in args.sequence_info
    ]

    #  Figure out the union of headers for all of these seq info files
    out_si_header = []
    for si_r in seq_info_readers:
        # Try to preserve order, only appending missing fn
        for fn in si_r.fieldnames:
            if fn not in out_si_header:
                out_si_header.append(fn)

    # We should at least have a seqname column in at least one seqinfo file.
    # If not, error and quit.
    if 'seqname' not in out_si_header:
        logging.error("No seqname column found in the seqinfo files. Cannot filter")
        sys.exit(-1)

    # Implicit else...
    out_h = args.output
    si_writer = csv.DictWriter(
        out_h,
        fieldnames=out_si_header)
    si_writer.writeheader()

    # Load all the seq_ids into a set
    seq_ids = set()
    for fasta_h in args.fasta:
        for sr in fastalite.fastalite(fasta_h):
            seq_ids.add(sr.id)

    logging.info("Found %d unique sequence IDs" % len(seq_ids))

    # Now go through each sequence information file. See if the row matches one of our target ids.
    for si_r in seq_info_readers:
        for row in si_r:
            if row['seqname'] in seq_ids:
                si_writer.writerow(row)
                seq_ids.remove(row['seqname'])
        # No need to continue if we have no seq ids to find
            if len(seq_ids) == 0:
                break
        if len(seq_ids) == 0:
                break

    # we should have found all of our seq IDs by now. Therefore express concern if we haven't

    if len(seq_ids) > 0:
        logging.error("Could not find sequence information for sequences in the FASTA file with IDs: %s" % (", ".join(seq_ids)))
        sys.exit(-1)

if __name__ == "__main__":
    main()
