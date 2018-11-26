#!/usr/bin/env python
import argparse
import fastalite
import logging
import csv


# UC Format for searching. TSV
# 0. Record type: H, or N. H= Hit. N= No hit
# 1. Ordinal number of the target sequence (based on input order, starting from zero). Set to ?*? for N.
# 2. Sequence length. Set to ?*? for N.
# 3. Percentage of similarity with the target sequence. Set to ?*? for N.
# 4. Match orientation + or -. . Set to ?.? for N.
# 5. Not used, always set to zero for H, or ?*? for N.
# 6. Not used, always set to zero for H, or ?*? for N.
# 7. Compact representation of the pairwise alignment using the CIGAR for- mat (Compact Idiosyncratic Gapped Alignment Report): M (match), D (deletion) and I (insertion). The equal sign ?=? indicates that the query is identical to the centroid sequence. Set to ?*? for N.
# 8. Label of the query sequence.
# 9. Label of the target centroid sequence. Set to ?*? for N.

def read_uc(uc_h):
    uc_reader = csv.DictReader(
        uc_h,
        delimiter='\t',
        fieldnames=[
            'result_code',
            'target_n',
            'seq_len',
            'percent_id',
            'match_orientation',
            'r5',
            'r6',
            'cigar',
            'query_id',
            'target_id',
        ])
    return [r for r in uc_reader]


def main():
    args_parser = argparse.ArgumentParser(
        description="""Filters sequences whose best search results (uc format) falls below a minimum
        percent sequence identity.
        """)

    args_parser.add_argument(
        'query_fasta',
        help='FASTA file containing all of the query sequences.',
        type=argparse.FileType(mode='r')
    )
    args_parser.add_argument(
        '--uc',
        nargs='+',
        help='UC files(s) with search results for these queries',
        type=argparse.FileType(mode='r'),
        required=True
    )
    args_parser.add_argument(
        '--min-best',
        '-m',
        help="""Minimum best identity (0.0 to 1.0)
        between query and reference seqs
        for a query to be considered matched""",
        required=True,
        type=float,
    )
    args_parser.add_argument(
        '--output',
        '-o',
        help="""FASTA file into which we should place
        our query sequences without a hit above minbest""",
        required=True,
        type=argparse.FileType('w')
    )

    args = args_parser.parse_args()
    logging.basicConfig(level=logging.INFO)

    min_best = float(args.min_best)
    out_h = args.output

    all_uc_data = []
    uc_h = args.uc

    for uc in uc_h:
        all_uc_data += read_uc(uc)

    logging.info("%d query result rows read in from the UC file(s)" % len(all_uc_data))

    # passed_queries: query sequence ids with a hit >= minbest
    passed_queries = set()

    # Get the queries
    query_ids = {r['query_id'] for r in all_uc_data}
    logging.info("%d unique query_ids searched." % len(query_ids))

    # For each query...
    for query_id in query_ids:
        # Grab rows where we have a hit
        query_hits = [
            r for r in all_uc_data
            if r['query_id'] == query_id if r['result_code'] == 'H'
        ]
        if len(query_hits) == 0:  # No hits for this query.
            continue
        if max([float(r['percent_id']) / 100.0 for r in query_hits]) >= min_best:
            passed_queries.add(query_id)

    logging.info("{} query_ids had a best hit meeting our threshold of {}.".format(
        len(passed_queries),
        min_best
    ))

    for sr in fastalite.fastalite(args.query_fasta):
        if sr.id in passed_queries:  # If we are in our passed queries, leave it out.
            continue
        # Implicit else
        if sr.id not in query_ids:
            logging.warn("%s was in the input query fasta but had no entry in the UC files. Included in the output" % sr.id)
        out_h.write(">%s %s\n%s\n" % (sr.id, sr.description, sr.seq))

if __name__ == "__main__":
    main()
