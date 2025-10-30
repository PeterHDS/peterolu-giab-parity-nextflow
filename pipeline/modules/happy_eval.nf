process HAPPY_EVAL {
  publishDir("${params.outdir}/happy", mode: 'copy')
  container "pkrusche/hap.py:latest"

  input:
    path refdir
    tuple path(truth_vcf_gz), path(regions_bed)
    tuple val(sample_id), path(call_vcf_gz)

  output:
    path "happy/*"

  script:
  """
  cp -r $refdir ./ref
  mkdir -p happy
  hap.py -r ref/ref.fa -f ${regions_bed} -o happy/${sample_id} ${truth_vcf_gz} ${call_vcf_gz}
  """

  stub:
  """
  mkdir -p happy
  echo "metric,value" > happy/${sample_id}.summary.csv
  echo "TP,0"        >> happy/${sample_id}.summary.csv
  echo "FP,0"        >> happy/${sample_id}.summary.csv
  echo "FN,0"        >> happy/${sample_id}.summary.csv
  """
}
