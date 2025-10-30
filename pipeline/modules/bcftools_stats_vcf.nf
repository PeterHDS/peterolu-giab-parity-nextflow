process BCFTOOLS_STATS_VCF {
  tag "$sample_id"
  publishDir("${params.outdir}/vcf_metrics", mode: 'copy')
  container "quay.io/biocontainers/bcftools:1.17--h3cc50cf_1"

  input:
    tuple val(sample_id), path(vcf_gz)

  output:
    path "${sample_id}.bcftools_stats.txt"

  script:
  """
  bcftools stats ${vcf_gz} > ${sample_id}.bcftools_stats.txt
  """
}
