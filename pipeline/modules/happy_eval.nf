process HAPPY_EVAL {
  publishDir("${params.outdir}/happy", mode: 'copy')
  container "mambaorg/micromamba:1.5.8"

  input:
    path refdir
    tuple path(truth_vcf_gz), path(regions_bed)
    tuple val(sample_id), path(call_vcf_gz)

  output:
    path "happy/*"

  shell:
  '''
  set -euo pipefail
  mkdir -p happy

  micromamba create -y -n happy -c conda-forge -c bioconda hap.py=0.3.15 >&2
  micromamba run -n happy hap.py \
    -r "!{refdir}/ref.fa" \
    -f "!{regions_bed}" \
    -o "happy/!{sample_id}" \
    "!{truth_vcf_gz}" \
    "!{call_vcf_gz}"
  '''
}
