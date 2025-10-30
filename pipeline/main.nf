nextflow.enable.dsl=2

include { FASTQC }         from './modules/fastqc.nf'
include { BWA_ALIGN }      from './modules/bwa_align.nf'
include { SAMTOOLS_SORT }  from './modules/samtools_sort.nf'
include { SAMTOOLS_STATS } from './modules/samtools_stats.nf'
include { MULTIQC }        from './modules/multiqc.nf'

// sample_id → [R1, R2]
read_pairs = Channel
  .fromFilePairs( params.reads )
  .map { sid, files -> tuple( sid, files ) }

ref_ch = Channel.value( file(params.reference_fasta) )

workflow {
  fastqc_out  = FASTQC(read_pairs)
  aligned_sam = BWA_ALIGN(ref_ch, read_pairs)
  sorted_bam  = SAMTOOLS_SORT(aligned_sam)
  stats_out   = SAMTOOLS_STATS(sorted_bam)
  MULTIQC( fastqc_out.mix(stats_out).collect() )
}
