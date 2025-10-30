nextflow.enable.dsl=2

include { FASTQC }  from './modules/fastqc.nf'
include { MULTIQC } from './modules/multiqc.nf'

// fromFilePairs WITHOUT flat:true so we get: [ id, [R1,R2] ]
Channel
  .fromFilePairs( params.reads )
  .map { sid, files -> tuple(sid, files) }
  .set { read_pairs }

workflow {
  // 1) Run FastQC -> emits ONE channel of directories ("fastqc")
  fastqc_dirs = FASTQC(read_pairs)

  // 2) Bundle all FastQC dirs into a single value for MultiQC
  fastqc_dirs.collect().set { fastqc_bundle }

  // 3) Aggregate with MultiQC
  MULTIQC(fastqc_bundle)
}
