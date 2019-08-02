#!/usr/bin/env python
import argparse
import logging
import fastalite

#
#   Given at least set(s) of paired reads in fastq format,
#   combine into one pair of reads, in order for each pair
#   Modest restrictions to keep things moving:
#       Pairs should be in order
#       Reads should be in order.
#


def get_seq_id(raw_id, normalize=True):
    if normalize:
        return raw_id.split('/')[0]
    else:
        return raw_id


def write_fastq(sr, out_h):
    out_h.write("@{} {}\n{}\n+\n{}\n".format(
        sr.id,
        sr.description,
        sr.seq,
        sr.qual
    ))


def main():
    args_parser = argparse.ArgumentParser(
        description="""Given set(s) of paired reads in fastq format
        combine all into one pair of reads also in fastq format.
        Concurrently confirm all R1 have matched R2.
        """
    )

    args_parser.add_argument(
        '--in-1',
        '-1',
        help='Read 1 Files to be combined.',
        nargs='+',
        required=True,
        type=fastalite.Opener(mode='r')
    )
    args_parser.add_argument(
        '--in-2',
        '-2',
        help="""Read 2 Files to be combined. Must be in same order as --in-1""",
        nargs='+',
        required=True,
        type=fastalite.Opener(mode='r')
    )

    args_parser.add_argument(
        '--out-1',
        '-o1',
        help='File into which we should place our combined R1',
        required=True,
        type=fastalite.Opener(mode='w')
    )
    args_parser.add_argument(
        '--out-2',
        '-o2',
        help='File into which we should place our combined R2',
        required=True,
        type=fastalite.Opener(mode='w')
    )

    args_parser.add_argument(
        '--normalize-ids',
        '-ni',
        help='Normalize IDs for pairs by stripping /x from the end',
        action='store_true'
    )

    args = args_parser.parse_args()
    logging.basicConfig(level=logging.INFO)

    assert len(args.in_1) == len(args.in_2), "Mismatched number of forward and reverse read files."

    # Loop 1: Identify ALL R1 and R2 IDs in all files.
    # Also look for duplicated IDs
    IDs_R1 = set()
    IDs_R2 = set()
    logging.info("Looping through files to identify all sequence IDs")
    for r1_h, r2_h in zip(args.in_1, args.in_2):
        # Will not use list comprehension to allow for error handling...
        file_ids_r1 = set()
        r1_reader = fastalite.fastqlite(r1_h)
        try:
            while True:
                try:
                    sr = next(r1_reader)
                    file_ids_r1.add(get_seq_id(sr.id, args.normalize_ids))
                except ValueError:
                    pass
        except StopIteration:
            pass

        file_ids_r2 = set()
        r2_reader = fastalite.fastqlite(r2_h)
        try:
            while True:
                try:
                    sr = next(r2_reader)
                    file_ids_r2.add(get_seq_id(sr.id, args.normalize_ids))
                except ValueError:
                    pass
        except StopIteration:
            pass

        if len(IDs_R1.intersection(file_ids_r1)) > 0:
            logging.warning("{:,} of {:,} R1 read IDs from this file overlap with others".format(
                len(IDs_R1.intersection(file_ids_r1)),
                len(file_ids_r1)
            ))
        if len(IDs_R2.intersection(file_ids_r2)) > 0:
            logging.warning("{:,} of {:,} R2 read IDs from this file overlap with others".format(
                len(IDs_R2.intersection(file_ids_r2)),
                len(file_ids_r2)
            ))
        IDs_R1.update(file_ids_r1)
        IDs_R2.update(file_ids_r2)
        r1_h.seek(0)
        r2_h.seek(0)

    overlapped_ids = IDs_R1.intersection(IDs_R2)
    starting_num_ids = len(overlapped_ids)
    logging.info("There are {:,} overlapping IDs from {:,} forward read IDs and {:,} reverse read IDs".format(
        starting_num_ids,
        len(IDs_R1),
        len(IDs_R2)
    ))

    for r1_h, r2_h in zip(args.in_1, args.in_2):
        srs_r1 = fastalite.fastqlite(r1_h)
        srs_r2 = fastalite.fastqlite(r2_h)
        sr_1 = None
        sr_2 = None
        try:
            while sr_1 is None:
                try:
                    sr_1 = next(srs_r1)
                except ValueError:
                    sr_1 = None
            while sr_2 is None:
                try:
                    sr_2 = next(srs_r2)
                except ValueError:
                    sr_2 = None
            while len(overlapped_ids) > 0:
                if len(overlapped_ids) % 10000 == 0:
                    logging.info(
                        "{:,} of {:,} pairs remaining".format(
                            len(overlapped_ids),
                            starting_num_ids,
                        )
                    )
                while (sr_1 is None) or (get_seq_id(sr_1.id, args.normalize_ids) not in overlapped_ids):
                    try:
                        sr_1 = next(srs_r1)
                    except ValueError:
                        sr_1 = None
                while (sr_2 is None) or (get_seq_id(sr_2.id, args.normalize_ids) not in overlapped_ids):
                    try:
                        sr_2 = next(srs_r2)
                    except ValueError:
                        sr_2 = None
                assert get_seq_id(sr_1.id, args.normalize_ids) == get_seq_id(sr_2.id, args.normalize_ids), "Order off of reads"
                # Implicit else paired and shared.
                # Remove it from the target list (takes care of duplicates)
                overlapped_ids.remove(get_seq_id(sr_1.id, args.normalize_ids))
                # Write out pair...
                write_fastq(sr_1, args.out_1)
                write_fastq(sr_2, args.out_2)
                # move to next
                try:
                    sr_1 = next(srs_r1)
                except ValueError:
                    pass
                try:
                    sr_2 = next(srs_r2)
                except ValueError:
                    pass
        except StopIteration:
            pass

if __name__ == "__main__":
    main()
