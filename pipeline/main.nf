nextflow.enable.dsl = 2

//
// ===== Parameters (defaults; overridden by nextflow.config if set) =====
//
params.reads            = params.reads            ?: 'data/*_{R1,R2}.fastq.gz'
params.reference_fasta  = params.reference_fasta  ?: 'ref/reference.fa'
params.truth_vcf        = params.truth_vcf        ?: 'truth/truth.vcf.gz'
params.confident_bed    = params.confident_bed    ?: 'truth/confident.bed'
params.outdir           = params.outdir           ?: './results'

//
// ===== Module imports / aliases =====
// Adjust module filenames if they differ in your repo
//

// Reference prep (faidx + dict)
include { REF_FAIDX as REF_FAIDX_A } from './modules/ref_faidx.nf'
include { REF_DICT  as REF_DICT_A  } from './modules/ref_dict.nf'

include { REF_FAIDX as REF_FAIDX_B } from './modules/ref_faidx.nf'
include { REF_DICT  as REF_DICT_B  } from './modules/ref_dict.nf'

// FastQC QC
include { FASTQC as FASTQC_A } from './modules/fastqc.nf'
include { FASTQC as FASTQC_B } from './modules/fastqc.nf'

// BWA alignment
include { BWA_ALIGN as BWA_ALIGN_A } from './modules/bwa_align.nf'
include { BWA_ALIGN as BWA_ALIGN_B } from './modules/bwa_align.nf'

// Samtools sort / stats / index
include { SAMTOOLS_SORT  as SAMTOOLS_SORT_A  } from './modules/samtools_sort.nf'
include { SAMTOOLS_SORT  as SAMTOOLS_SORT_B  } from './modules/samtools_sort.nf'

include { SAMTOOLS_STATS as SAMTOOLS_STATS_A } from './modules/samtools_stats.nf'
include { SAMTOOLS_STATS as SAMTOOLS_STATS_B } from './modules/samtools_stats.nf'

include { SAMTOOLS_INDEX as SAMTOOLS_INDEX_A } from './modules/samtools_index.nf'
include { SAMTOOLS_INDEX as SAMTOOLS_INDEX_B } from './modules/samtools_index.nf'

// GATK HaplotypeCaller
include { GATK_HC as GATK_HC_A } from './modules/gatk_hc.nf'
include { GATK_HC as GATK_HC_B } from './modules/gatk_hc.nf'

// bcftools VCF stats (Branch A only for now)
include { BCFTOOLS_STATS_VCF } from './modules/bcftools_stats_vcf.nf'

// MultiQC (no A/B split, just a global report)
include { MULTIQC } from './modules/multiqc.nf'

//
// ========================= Main workflow ===============================
//
workflow {

  //
  // Inputs
  //
  read_pairs_A = Channel
    .fromFilePairs( params.reads )
    .map { sid, files -> tuple(sid, files) }

  read_pairs_B = Channel
    .fromFilePairs( params.reads )
    .map { sid, files -> tuple(sid, files) }

  ref_in_A = Channel.value( file(params.reference_fasta) )
  ref_in_B = Channel.value( file(params.reference_fasta) )

  // GIAB truth kept defined, but unused for now
  truth_chan = Channel.value(
    tuple(
      file(params.truth_vcf),
      file(params.confident_bed)
    )
  )

  //
  // ===== Branch A: QC → map → sort → stats → index → HC → bcftools stats =====
  //
  refA        = REF_FAIDX_A(ref_in_A)
  ref_ready_A = REF_DICT_A(refA)    // emits prepared reference

  fastqc_A    = FASTQC_A(read_pairs_A)
  aligned_A   = BWA_ALIGN_A(ref_in_A, read_pairs_A)
  sorted_A    = SAMTOOLS_SORT_A(aligned_A)
  stats_A     = SAMTOOLS_STATS_A(sorted_A)
  indexed_A   = SAMTOOLS_INDEX_A(sorted_A)   // -> (sid, bam, bai)

  vcf_A_raw   = GATK_HC_A(ref_ready_A, indexed_A)   // -> (sid, vcf.gz)

  // extract just the VCF path for bcftools stats
  vcf_A_only  = vcf_A_raw.map { sid, vcf -> vcf }
  vcf_stats_A = BCFTOOLS_STATS_VCF( vcf_A_only )

  //
  // ===== Branch B: QC → map → sort → stats → index → HC =====
  // (hap.py evaluation temporarily disabled)
  //
  refB        = REF_FAIDX_B(ref_in_B)
  ref_ready_B = REF_DICT_B(refB)

  fastqc_B    = FASTQC_B(read_pairs_B)
  aligned_B   = BWA_ALIGN_B(ref_in_B, read_pairs_B)
  sorted_B    = SAMTOOLS_SORT_B(aligned_B)
  stats_B     = SAMTOOLS_STATS_B(sorted_B)
  indexed_B   = SAMTOOLS_INDEX_B(sorted_B)    // -> (sid, bam, bai)

  vcf_B_raw   = GATK_HC_B(ref_ready_B, indexed_B)   // -> (sid, vcf.gz)

  // NOTE: no HAPPY_EVAL_PROC() call here for now

  //
  // MultiQC – use Branch A outputs as demo inputs
  //
  multiqc_input = fastqc_A
                    .mix(stats_A)
                    .mix(vcf_stats_A)

  MULTIQC( multiqc_input )
}
