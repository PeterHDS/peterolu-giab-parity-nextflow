process REF_DICT {
  publishDir("${params.outdir}/ref", mode: 'copy')
  container "quay.io/biocontainers/gatk4:4.4.0.0--py39hdfd78af_0"

  input:
    path refdir

  output:
    path "ref"

  script:
  """
  cp -r $refdir ./ref
  gatk CreateSequenceDictionary -R ref/ref.fa -O ref/ref.dict
  """

  stub:
  """
  mkdir -p ref
  echo "stub dict" > ref/ref.dict
  """
}
