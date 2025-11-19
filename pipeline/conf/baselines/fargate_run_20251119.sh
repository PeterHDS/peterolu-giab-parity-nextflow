#!/usr/bin/env bash
set -euo pipefail

conda activate giab-nf

export AWS_REGION=eu-west-2
export NXF_S3=1
export BUCKET=s3://portable-precision-peterhds-eu2
export QUEUE_FG=pp-queue-fg
export BATCH_JOB_ROLE_ARN=arn:aws:iam::556008108909:role/pp-batch-task-role
export BATCH_EXEC_ROLE_ARN=arn:aws:iam::556008108909:role/pp-batch-exec-role

nextflow run main.nf \
  -profile awsbatch_fg \
  -work-dir "${BUCKET}/work" \
  --outdir "${BUCKET}/results/fargate-$(date +%Y%m%d-%H%M%S)" \
  -with-report awsbatch_report.html \
  -with-timeline awsbatch_timeline.html \
  -with-trace awsbatch_trace.txt
