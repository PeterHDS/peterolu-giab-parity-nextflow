nextflow.enable.dsl=2

include { FASTQC } from './modules/fastqc.nf'

// Build tuples: [ sample_id, R1.fastq.gz, R2.fastq.gz ]
Channel
  .fromFilePairs( params.reads, size: 2 )         // returns [id, [R1, R2]]
  .map { id, pair -> tuple(id, pair[0], pair[1]) } // make it a 3-tuple
  .set { read_pairs }

workflow {
  FASTQC(read_pairs)
}
