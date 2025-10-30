process SAMTOOLS_SORT {
  tag "$sample_id"
  publishDir("${params.outdir}/bam", mode: 'copy')
  container "quay.io/biocontainers/samtools:1.17--hd87286a_1"

  input:
    tuple val(sample_id), path(sam)

  output:
    tuple val(sample_id), path("${sample_id}.sorted.bam")

  script:
  """
  samtools view -bS ${sam} | samtools sort -o ${sample_id}.sorted.bam
  samtools index ${sample_id}.sorted.bam
  """

  stub:
  """
  # minimal BAM placeholder
  : > ${sample_id}.sorted.bam
  """
}
