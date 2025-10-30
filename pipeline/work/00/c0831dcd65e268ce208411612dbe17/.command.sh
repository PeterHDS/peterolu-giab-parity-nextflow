#!/bin/bash -ue
samtools view -bS demo.sam | samtools sort -o demo.sorted.bam
samtools index demo.sorted.bam
