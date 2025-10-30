process MULTIQC {
  publishDir("${params.outdir}/multiqc", mode: 'copy')
  container "quay.io/biocontainers/multiqc:1.19--pyhdfd78af_0"

  input:
    path inputs

  output:
    path "multiqc/*"

  script:
  """
  mkdir -p multiqc
  multiqc -o multiqc .
  """

  stub:
  """
  mkdir -p multiqc
  echo "stub multiqc output" > multiqc/summary.txt
  """
}
