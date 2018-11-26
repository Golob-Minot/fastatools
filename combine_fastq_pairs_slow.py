#!/usr/bin/env python
import argparse
import re
import logging
import fastalite
from collections import namedtuple

#
#   Given at least set(s) of paired reads in fastq format,
#   combine into one pair of reads, in order for each pair
#   this is an ultra-slow implementation that is super robust
#

re_preSeq = re.compile(r'^@(?P<id>\S+)( (?P<description>.+)|())')
re_preQual = re.compile(r'^\+((?P<id>\S+)( (?P<description>.+)|())|())')


def get_seq_id(raw_id, normalize=True):
    if normalize:
        return raw_id.split('/')[0]
    else:
        return raw_id


def parse_fastq_lines(lines, validate_annotations=False):
    #  Given a block of lines corresponding to a fastq record
    #  parse out four variables:
    #   - id
    #   - description
    #   - seq
    #   - qual

    rec_id = None
    rec_description = None
    rec_seq = ""
    rec_qual = ""

    rec_m = re_preSeq.match(lines[0])
    assert rec_m is not None, "Record does not start with an @id description line, as per fastq format"
    rec_groupdict = rec_m.groupdict()
    rec_id = rec_groupdict.get('id')
    rec_description = rec_groupdict.get('description')

    i = 1
    while i < len(lines) and lines[i][0] != '+':
        rec_seq += lines[i]
        i += 1
    # i is now at the qual identifier
    if validate_annotations:
        qual_intro_groupdict = re_preQual.match(lines[i]).groupdict()
        if qual_intro_groupdict.get('id', None) is not None:
            assert rec_groupdict.get('id') == qual_intro_groupdict.get('id')
        if qual_intro_groupdict.get('description', None) is not None:
            assert rec_groupdict.get('description') == qual_intro_groupdict.get('description')
    # Jump to where qual starts
    i += 1
    while i < len(lines):
        rec_qual += lines[i]
        i += 1
    record = namedtuple(
        "fastq_record",
        ['id', 'description', 'seq', 'qual']
    )
    record.id = rec_id
    if rec_description is not None:
        record.description = rec_description
    else:
        record.description = ""
    record.seq = rec_seq
    record.qual = rec_qual
    return record


def main():
    args_parser = argparse.ArgumentParser(
        description="""Given set(s) of paired reads in fastq format
        combine all into one pair of reads also in fastq format.
        Concurrently confirm all R1 have matched R2
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
        help="""Read 2 Files to be combined.""",
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

    # A mapping of seq IDs to a tuple (file_h, tell_position)
    seq_ids_to_r1 = {}
    seq_ids_to_r2 = {}

    # Loop through forward reads to get all the IDs
    logging.info(
        "Reading in sequence IDs from %d forward read files",
        len(args.in_1)
    )
    #
    #  FASTQ format is in four line blocks:
    #   0. @id description
    #   1. sequence
    #   2. +id description (optional beyond +)
    #   3. qual
    #
    #   Sadly '@' is a valid character in qual lines (UGGGH)
    #   so we cannot just search for that.
    #

    for r1_h in args.in_1:
        # Make an mmap of the file
        r1_mm = mmap.mmap(r1_h., 0, prot=mmap.PROT_READ)
        rec_block = None
        cur_pos = r1_mm.tell()
        line = r1_mm.readline().decode('utf-8')
        print(line)
        return
        rec_line_n = 0
        while line:
            if rec_line_n == 0:  # New record
                # Handle previous record here
                if rec_block is not None:
                    record = parse_fastq_lines(rec_block)
                    seq_ids_to_r1[get_seq_id(record.id, args.normalize_ids)] = (
                        r1_h,
                        block_start
                    )
                # Now create a new record block
                rec_block = [line.strip()]
                block_start = cur_pos
            else:
                rec_block.append(line.strip())
                cur_pos = r1_h.tell()
            rec_line_n += 1
            if rec_line_n == 4:
                rec_line_n = 0
            line = r1_h.readline()

    logging.info(
        "Reading in sequence IDs from %d reverse read files",
        len(args.in_2)
    )
    for r2_h in args.in_1:
        r2_h.seek(0)
        rec_block = None
        cur_pos = r2_h.tell()
        line = r2_h.readline()
        rec_line_n = 0
        while line:
            if rec_line_n == 0:  # New record
                # Handle previous record here
                if rec_block is not None:
                    record = parse_fastq_lines(rec_block)
                    seq_ids_to_r2[get_seq_id(record.id, args.normalize_ids)] = (
                        r2_h,
                        block_start
                    )
                # Now create a new record block
                rec_block = [line.strip()]
                block_start = cur_pos
            else:
                rec_block.append(line.strip())
                cur_pos = r2_h.tell()
            rec_line_n += 1
            if rec_line_n == 4:
                rec_line_n = 0
            line = r2_h.readline()

    # Use sets to get the IDs shared between all the forward and reverse reads
    shared_ids = list(set(seq_ids_to_r1.keys()).intersection(set(seq_ids_to_r2.keys())))
    logging.info(
        "{:,} shared IDs from {:,} R1 IDs and {:,} R2 IDs".format(
            len(shared_ids),
            len(seq_ids_to_r1),
            len(seq_ids_to_r2)
        )
    )
    logging.info("Writing R1 and R2 to combined fastq file")
    for n, seq_id in enumerate(shared_ids):
        if (n + 1) % 100 == 0:
            logging.info("Working on seq id {:,} of {:,}".format(n + 1, len(shared_ids)))
        r1_h, r1_pos = seq_ids_to_r1[seq_id]
        r1_h.seek(r1_pos)
        rec_block_r1 = [
            r1_h.readline().strip(),
            r1_h.readline().strip(),
            r1_h.readline().strip(),
            r1_h.readline().strip()
        ]
        sr_r1 = parse_fastq_lines(rec_block_r1)

        r2_h, r2_pos = seq_ids_to_r2[seq_id]
        r2_h.seek(r2_pos)
        rec_block_r2 = [
            r2_h.readline().strip(),
            r2_h.readline().strip(),
            r2_h.readline().strip(),
            r2_h.readline().strip()
        ]
        sr_r2 = parse_fastq_lines(rec_block_r2)

        try:
            assert get_seq_id(sr_r1.id, args.normalize_ids) == get_seq_id(sr_r2.id, args.normalize_ids), "IDs mismatch in file for {}".format(seq_id)
        except AssertionError as e:
            logging.error("{}".format(e))
            continue
        if sr_r1.seq == "" or sr_r2.seq == "":
            logging.warning("SeqID {} has an empty sequence.".format(seq_id))
            continue

        args.out_1.write("@{} {}\n{}\n+\n{}\n".format(
            sr_r1.id,
            sr_r1.description,
            sr_r1.seq,
            sr_r1.qual
        ))
        args.out_2.write("@{} {}\n{}\n+\n{}\n".format(
            sr_r2.id,
            sr_r2.description,
            sr_r2.seq,
            sr_r2.qual
        ))


if __name__ == "__main__":
    main()
