# fastatools
#
# VERSION               golob/fastatools:0.7.1__bcw.0.3.1

FROM      alpine:3.8
RUN apk add --no-cache python3-dev build-base bash
RUN ln -s /usr/bin/python3 /usr/bin/python
RUN pip3 install \
boto3 \
awscli \
biopython \
bucket_command_wrapper==0.3.1 \
fastalite>=0.3

ADD combine_fasta.py /usr/local/bin
ADD combine_fastq_pairs.py /usr/local/bin
ADD fasta_a_not_b.py /usr/local/bin
ADD fasta_seq_info.py /usr/local/bin
ADD seqs_below_minbest.py /usr/local/bin

RUN chmod +x /usr/local/bin/*.py

RUN mkdir /fh && mkdir /app
RUN mkdir -p /mnt/inputs/file && mkdir -p /mnt/outputs/file && mkdir /scratch
