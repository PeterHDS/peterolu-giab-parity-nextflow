process BCFTOOLS_STATS_VCF {
  tag "$sample_id"
  publishDir("${params.outdir}/vcf_metrics", mode: 'copy')
  container "quay.io/biocontainers/bcftools:1.17--hf3cf87c_0"

  input:
    tuple val(sample_id), path(vcf_gz)

  output:
    path "${sample_id}.bcftools.stats.txt"

  script:
  """
  bcftools stats ${vcf_gz} > ${sample_id}.bcftools.stats.txt
  """

  stub:
  """
  echo "bcftools stats (stub)" > ${sample_id}.bcftools.stats.txt
  """
}
