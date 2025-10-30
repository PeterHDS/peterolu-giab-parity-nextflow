process FASTQC {
  tag "$sample_id"
  publishDir("${params.outdir}/fastqc", mode: 'copy')
  container "quay.io/biocontainers/fastqc:0.12.1--hdfd78af_0"

  input:
    tuple val(sample_id), path(reads)

  output:
    path "fastqc"

  script:
  """
  mkdir -p fastqc
  fastqc -o fastqc $reads
  """

  stub:
  """
  mkdir -p fastqc
  echo "stub fastqc output" > fastqc/${sample_id}_fastqc.txt
  """
}
