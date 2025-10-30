#!/bin/bash -ue
mkdir -p fastqc
fastqc -o fastqc demo_R1.fastq.gz demo_R2.fastq.gz
