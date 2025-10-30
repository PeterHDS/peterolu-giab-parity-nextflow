process REF_DICT {
  container "quay.io/biocontainers/gatk4:4.2.5.0--hdfd78af_0"

  input:
    path refdir

  output:
    path "ref"

  script:
  """
  cp -r $refdir ./ref
  gatk CreateSequenceDictionary -R ref/ref.fa -O ref/ref.dict
  """
}
