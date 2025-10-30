process SAMTOOLS_STATS {
  tag "$sample_id"
  publishDir("${params.outdir}/metrics", mode: 'copy')
  container "quay.io/biocontainers/samtools:1.17--hd87286a_1"

  input:
    tuple val(sample_id), path(bam)

  output:
    path "${sample_id}.samtools_flagstat.txt"
    path "${sample_id}.samtools_stats.txt"

  script:
  """
  samtools flagstat ${bam} > ${sample_id}.samtools_flagstat.txt
  samtools stats ${bam} > ${sample_id}.samtools_stats.txt
  """

  stub:
  """
  echo "This file was produced by samtools flagstat (stub)" > ${sample_id}.samtools_flagstat.txt
  echo "This file was produced by samtools stats (stub)"   > ${sample_id}.samtools_stats.txt
  """
}
