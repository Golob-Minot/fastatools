# fastatools
#
# VERSION               golob/fastatools:0.8.0A

FROM quay.io/biocontainers/biopython:1.78

RUN ln -s /usr/bin/python3 /usr/bin/python

RUN pip3 install \
fastalite>=0.3

ADD combine_fasta.py /usr/local/bin
ADD combine_fastq_pairs.py /usr/local/bin
ADD fasta_a_not_b.py /usr/local/bin
ADD fasta_seq_info.py /usr/local/bin
ADD seqs_below_minbest.py /usr/local/bin

RUN chmod +x /usr/local/bin/*.py