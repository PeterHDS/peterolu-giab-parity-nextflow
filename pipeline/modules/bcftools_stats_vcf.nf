process BCFTOOLS_STATS_VCF {
  tag "${vcf.simpleName}"
  publishDir "${params.outdir}/vcf_stats", mode: 'copy'
  container 'quay.io/biocontainers/bcftools:1.17--h3cc50cf_1'

  input:
    // single VCF(.gz) path (no tuples!)
    path vcf

  output:
    path "${vcf.simpleName}.bcftools.stats.txt", emit: stats

  script:
  """
  set -euo pipefail
  # index if missing (harmless if already present)
  if [ ! -f "${vcf}.tbi" ] && [ ! -f "${vcf}.csi" ]; then
    bcftools index -f "${vcf}"
  fi

  bcftools stats "${vcf}" > "${vcf.simpleName}.bcftools.stats.txt"
  """
}
