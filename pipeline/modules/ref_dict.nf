process REF_DICT {
  tag "Create dict"
  publishDir("${params.outdir}/ref", mode: 'copy')
  container 'broadinstitute/gatk:4.4.0.0'

  input:
    // Stage the incoming reference directory AS 'ref' so the path is stable and no self-copies happen
    path refdir, stageAs: 'ref'

  output:
    // Re-emit the prepared dir for downstream steps
    path 'ref'

  script:
  """
  set -euo pipefail
  gatk CreateSequenceDictionary -R ref/ref.fa -O ref/ref.dict
  """
}
