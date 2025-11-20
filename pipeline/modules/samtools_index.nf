process SAMTOOLS_INDEX {
  tag "$sid"
  container 'staphb/samtools:1.19'

  input:
    tuple val(sid), path(bam)

  output:
    // emit bam and its bai together so downstream tasks receive both
    tuple val(sid), path(bam), path("${bam}.bai")

  script:
  """
  set -euo pipefail
  samtools index -b ${bam}
  """
}
