process FASTQC {
  tag "$sample_id"
  publishDir("${params.outdir}/fastqc", mode: 'copy')
  container 'quay.io/biocontainers/fastqc:0.12.1--hdfd78af_0'

  input:
    tuple val(sample_id), path(read1), path(read2)

  output:
    path "fastqc/*_fastqc.html"
    path "fastqc/*_fastqc.zip"

  script:
  """
  mkdir -p fastqc
  fastqc -o fastqc $read1 $read2
  """
}
