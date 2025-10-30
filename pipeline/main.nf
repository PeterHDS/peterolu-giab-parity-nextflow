nextflow.enable.dsl=2

process HELLO {
  publishDir(params.outdir), mode: 'copy'
  """
  echo 'hello portable precision' > hello.txt
  """
}

workflow {
  HELLO()
}
