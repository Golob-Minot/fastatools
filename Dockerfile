# fastatools
#
# VERSION               golob/fastatools:0.7.1__bcw.0.3.1A

FROM biocontainers/biocontainers:v1.1.0_cv2
USER root
RUN apt update && apt install -y python3 python3-pip
RUN ln -s /usr/bin/python3 /usr/bin/python
USER biodocker
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

USER root
RUN chmod +x /usr/local/bin/*.py

RUN mkdir /fh && mkdir /app
RUN mkdir -p /mnt/inputs/file && mkdir -p /mnt/outputs/file && mkdir /scratch
USER biodocker
