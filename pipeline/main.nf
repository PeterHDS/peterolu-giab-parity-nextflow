nextflow.enable.dsl=2

include { FASTQC }             from './modules/fastqc.nf'
include { BWA_ALIGN }          from './modules/bwa_align.nf'
include { SAMTOOLS_SORT }      from './modules/samtools_sort.nf'
include { SAMTOOLS_STATS }     from './modules/samtools_stats.nf'
include { REF_FAIDX }          from './modules/ref_faidx.nf'
include { REF_DICT }           from './modules/ref_dict.nf'
include { GATK_HC }            from './modules/gatk_hc.nf'
include { BCFTOOLS_STATS_VCF } from './modules/bcftools_stats_vcf.nf'
include { MULTIQC }            from './modules/multiqc.nf'

// Pairs: [ sample_id, [R1,R2] ]
read_pairs = Channel
  .fromFilePairs( params.reads )
  .map { sid, files -> tuple( sid, files ) }

ref_in = Channel.value( file(params.reference_fasta) )

workflow {
  // Reference indices
  ref_faidx = REF_FAIDX(ref_in)
  ref_ready = REF_DICT(ref_faidx)

  // Read QC
  fastqc_out = FASTQC(read_pairs)

  // Align -> sort/index -> stats
  aligned_sam = BWA_ALIGN(ref_in, read_pairs)
  sorted_bam  = SAMTOOLS_SORT(aligned_sam)
  stats_out   = SAMTOOLS_STATS(sorted_bam)

  // Variant calling + VCF stats
  vcf_out   = GATK_HC(ref_ready, sorted_bam)
  vcf_stats = BCFTOOLS_STATS_VCF(vcf_out)

  // Aggregate
  MULTIQC( fastqc_out.mix(stats_out).mix(vcf_stats).collect() )
}
