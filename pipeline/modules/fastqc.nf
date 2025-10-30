process FASTQC {
  tag "$sample_id"
  publishDir("${params.outdir}/fastqc", mode: 'copy')
  container "quay.io/biocontainers/fastqc:0.12.1--hdfd78af_0"

  input:
    tuple val(sample_id), path(reads)   // reads is a LIST [R1,R2]

  output:
    path "fastqc"   // single directory per sample

  script:
  """
  mkdir -p fastqc
  fastqc -o fastqc $reads
  """
}
