process MULTIQC {

  tag "multiqc"

  // Publish MultiQC report into the main outdir
  publishDir "${params.outdir}/multiqc", mode: 'copy'

  input:
    // We don’t actually need to *use* qc_files in the script,
    // but this keeps the dependency so MultiQC runs at the end.
    path qc_files

  output:
    path "multiqc_report.html"
    path "multiqc_data", optional: true

  // Do NOT let MultiQC failure kill the whole pipeline
  errorStrategy 'ignore'
  maxRetries 1

  script:
  """
  mkdir -p multiqc

  # Run MultiQC but skip the problematic bcftools module
  multiqc --exclude bcftools -o multiqc .

  # Normalise outputs for Nextflow
  if [ -f multiqc/multiqc_report.html ]; then
    cp multiqc/multiqc_report.html .
  fi

  if [ -d multiqc/multiqc_data ]; then
    cp -r multiqc/multiqc_data .
  fi
  """
}
