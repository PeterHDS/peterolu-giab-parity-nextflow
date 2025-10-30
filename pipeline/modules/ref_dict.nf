process REF_DICT {
  container "quay.io/biocontainers/gatk4:4.2.5.0--hdfd78af_0"

  input:
    path refdir

  output:
    path "ref"

  shell:
  '''
  set -euo pipefail
  cp -r ref ./ref

  if [ -s ref/ref.dict ]; then
    echo "[REF_DICT] ref.dict already present; skipping." >&2
  else
    gatk CreateSequenceDictionary -R ref/ref.fa -O ref/ref.dict
  fi
  '''
}
