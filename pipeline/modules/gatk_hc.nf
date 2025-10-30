process GATK_HC {
  tag "$sample_id"
  publishDir("${params.outdir}/vcf", mode: 'copy')
  container "quay.io/biocontainers/gatk4:4.2.5.0--hdfd78af_0"

  input:
    path refdir
    tuple val(sample_id), path(bam)

  output:
    tuple val(sample_id), path("${sample_id}.vcf.gz")

  script:
  """
  cp -r $refdir ./ref
  gatk HaplotypeCaller \
    -R ref/ref.fa \
    -I ${bam} \
    -O ${sample_id}.vcf.gz
  """
}
