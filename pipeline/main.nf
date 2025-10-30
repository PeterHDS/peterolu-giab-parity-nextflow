nextflow.enable.dsl=2

// ----- Aliased includes so each branch can use its own copy -----
include { FASTQC          as FASTQC_A          } from './modules/fastqc.nf'
include { FASTQC          as FASTQC_B          } from './modules/fastqc.nf'
include { BWA_ALIGN       as BWA_ALIGN_A       } from './modules/bwa_align.nf'
include { BWA_ALIGN       as BWA_ALIGN_B       } from './modules/bwa_align.nf'
include { SAMTOOLS_SORT   as SAMTOOLS_SORT_A   } from './modules/samtools_sort.nf'
include { SAMTOOLS_SORT   as SAMTOOLS_SORT_B   } from './modules/samtools_sort.nf'
include { SAMTOOLS_STATS  as SAMTOOLS_STATS_A  } from './modules/samtools_stats.nf'
include { SAMTOOLS_STATS  as SAMTOOLS_STATS_B  } from './modules/samtools_stats.nf'
include { REF_FAIDX       as REF_FAIDX_A       } from './modules/ref_faidx.nf'
include { REF_FAIDX       as REF_FAIDX_B       } from './modules/ref_faidx.nf'
include { REF_DICT        as REF_DICT_A        } from './modules/ref_dict.nf'
include { REF_DICT        as REF_DICT_B        } from './modules/ref_dict.nf'
include { GATK_HC         as GATK_HC_A         } from './modules/gatk_hc.nf'
include { GATK_HC         as GATK_HC_B         } from './modules/gatk_hc.nf'
include { BCFTOOLS_STATS_VCF                 } from './modules/bcftools_stats_vcf.nf'
include { HAPPY_EVAL                          } from './modules/happy_eval.nf'
include { MULTIQC                             } from './modules/multiqc.nf'

// -------- Inputs (duplicate for two independent branches) --------
read_pairs_A = Channel
  .fromFilePairs( params.reads )
  .map { sid, files -> tuple( sid, files ) }

read_pairs_B = Channel
  .fromFilePairs( params.reads )
  .map { sid, files -> tuple( sid, files ) }

ref_in_A = Channel.value( file(params.reference_fasta) )
ref_in_B = Channel.value( file(params.reference_fasta) )

// GIAB toy truth bundle
truth_bundle = Channel.of( tuple(
  file(params.truth_vcf),
  file(params.confident_bed)
))

workflow {

  // ---- Branch A (stats path) ----
  ref_faidx_A   = REF_FAIDX_A(ref_in_A)
  ref_ready_A   = REF_DICT_A(ref_faidx_A)

  fastqc_A      = FASTQC_A(read_pairs_A)
  aligned_A     = BWA_ALIGN_A(ref_in_A, read_pairs_A)
  sorted_A      = SAMTOOLS_SORT_A(aligned_A)
  stats_map_A   = SAMTOOLS_STATS_A(sorted_A)

  vcf_A         = GATK_HC_A(ref_ready_A, sorted_A)
  vcf_stats_A   = BCFTOOLS_STATS_VCF(vcf_A)

  // ---- Branch B (hap.py path) ----
  ref_faidx_B   = REF_FAIDX_B(ref_in_B)
  ref_ready_B   = REF_DICT_B(ref_faidx_B)

  fastqc_B      = FASTQC_B(read_pairs_B)
  aligned_B     = BWA_ALIGN_B(ref_in_B, read_pairs_B)
  sorted_B      = SAMTOOLS_SORT_B(aligned_B)
  stats_map_B   = SAMTOOLS_STATS_B(sorted_B)

  vcf_B         = GATK_HC_B(ref_ready_B, sorted_B)
  happy_out     = HAPPY_EVAL(ref_ready_B, truth_bundle, vcf_B)

  // ---- Aggregate QC + mapping + VCF stats from Branch A only
  MULTIQC( fastqc_A.mix(stats_map_A).mix(vcf_stats_A).collect() )
}
