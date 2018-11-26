# fastatools
#
# VERSION               golob/fastatools:0.5__bcw.0.3.0

FROM      alpine:3.7
RUN mkdir /fh && mkdir /app
RUN mkdir -p /mnt/inputs/file && mkdir -p /mnt/outputs/file && mkdir /scratch
RUN mkdir -p /raxml/working
RUN apk add --no-cache python3 bash
RUN ln -s /usr/bin/python3 /usr/bin/python
RUN pip3 install \
boto3 \
awscli \
bucket_command_wrapper==0.3.0 \
fastalite>=0.3

ADD combine_fasta.py /usr/local/bin
ADD fasta_a_not_b.py /usr/local/bin
ADD fasta_seq_info.py /usr/local/bin
ADD seqs_below_minbest.py /usr/local/bin

RUN chmod +x /usr/local/bin/*.py