process REF_FAIDX {
  container "quay.io/biocontainers/samtools:1.17--hd87286a_1"

  input:
    path reference_fa

  output:
    path "ref"

  script:
  """
  mkdir -p ref
  cp $reference_fa ref/ref.fa
  samtools faidx ref/ref.fa
  """
}
