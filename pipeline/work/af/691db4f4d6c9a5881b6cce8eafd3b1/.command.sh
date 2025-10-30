#!/bin/bash -ue
samtools flagstat demo.sorted.bam > demo.flagstat.txt
samtools stats demo.sorted.bam > demo.stats.txt
