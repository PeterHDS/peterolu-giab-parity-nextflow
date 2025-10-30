#!/bin/bash -ue
cp miniref.fa ref.fa
bwa index ref.fa
bwa mem -t 2 ref.fa demo_R1.fastq.gz demo_R2.fastq.gz > demo.sam
