process GATK_HC {
  tag "$sid"
  // If you pinned resources in nextflow.config, those will override these
  cpus 2
  memory '6 GB'
  container 'broadinstitute/gatk:4.4.0.0'

  input:
    // stage the reference directory exactly as "ref" so path is stable
    path refdir, stageAs: 'ref'
    // CRITICAL: bring both BAM and BAI so NF stages both
    tuple val(sid), path(bam), path(bai)

  output:
    // emit sid + VCF path (bgzipped)
    tuple val(sid), path("${sid}.vcf.gz")

  script:
  """
  set -euo pipefail
  gatk HaplotypeCaller \
    -R ref/ref.fa \
    -I ${bam} \
    --native-pair-hmm-threads ${task.cpus ?: 2} \
    -O ${sid}.vcf.gz
  """
}
