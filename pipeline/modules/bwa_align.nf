process BWA_ALIGN {
  tag "$sample_id"
  publishDir("${params.outdir}/align", mode: 'copy')
  container "quay.io/biocontainers/bwa:0.7.17--he4a0461_11"

  input:
    path reference_fa
    tuple val(sample_id), path(reads)

  output:
    tuple val(sample_id), path("${sample_id}.sam")

  script:
  """
  cp $reference_fa ref.fa
  bwa index ref.fa
  bwa mem -t 2 -R "@RG\\tID:${sample_id}\\tSM:${sample_id}\\tPL:ILLUMINA" ref.fa $reads > ${sample_id}.sam
  """

  stub:
  """
  echo "@SQ\\tSN:stub\\tLN:100" > ${sample_id}.sam
  echo "@RG\\tID:${sample_id}\\tSM:${sample_id}" >> ${sample_id}.sam
  echo "@PG\\tID:stub" >> ${sample_id}.sam
  """
}
