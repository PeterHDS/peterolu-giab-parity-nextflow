# Portable Precision — Troubleshooting Log (by Peter Olu DS)

This is a short, honest log of how I took my tiny Nextflow pipeline from “nearly there” to “green ticks everywhere”. It’s written for humans, not robots.

## What I’m building
A small but production-style germline pipeline in Nextflow DSL2:
FastQC → BWA-MEM2 → samtools (sort/stats) → GATK HaplotypeCaller → MultiQC,
plus a hap.py check on a mock GIAB slice. Goal: show reproducibility and portability.

---

## What went wrong and how I fixed it

### 1) Job stuck at “0 of 1” (HAPPY_EVAL starved)
- Symptom: `HAPPY_EVAL (1) | 0 of 1` and nothing started.
- Why it happened: In DSL2, a channel feeds only one consumer. I was trying to feed the same data into two branches.
- Fix: I created two simple branches with their own inputs and **aliased** modules (e.g. `REF_DICT_A` / `REF_DICT_B`, `GATK_HC_A` / `GATK_HC_B`). No fancy operators needed.

### 2) “Process has been already used”
- Symptom: Nextflow complained that I reused the same module.
- Why: In DSL2, you can’t call the same included process twice unless you alias it.
- Fix: Aliased every module used in both branches.

### 3) Reference dictionary clashes + publish collisions
- Symptoms:
  - Picard/GATK: `ref/ref.dict already exists`
  - Nextflow: “Failed to publish … ref/ref”
- Why:
  - I generated the same dictionary twice and tried to publish to the same place twice.
- Fix:
  - Made `REF_DICT` **idempotent** (skip if file exists).
  - Stopped publishing this intermediate. Downstream steps read from work dirs anyway.

### 4) hap.py failed: output folder + VCF header
- Symptoms:
  - “Output path does not exist. Use -o …”
  - VCF header errors (missing `##FORMAT=GT`, missing `##contig`, etc.)
- Why:
  - I didn’t create the output folder before running hap.py.
  - My truth VCF header wasn’t valid for comparison.
- Fix:
  - `mkdir -p happy` inside the process before hap.py.
  - Rebuilt a tiny, valid truth VCF with proper headers, then bgzip + tabix indexed it.
  - BED region covers the variant; FASTA contig name matches.

### 5) Small but real: container tags + line endings
- I pinned known-good container tags and sanity-checked `--version`.
- I enforced LF line endings to stop weird Groovy parsing errors.

---

## Result
Everything ran cleanly:
- `HAPPY_EVAL` ✔ (hap.py produced summary CSVs)
- MultiQC ✔ (FastQC + samtools + VCF stats)
- End-to-end on mock data in ~2 minutes on my machine

---

## What I learned (in one breath)
Design branches on purpose, alias when you reuse modules, make reference/index steps idempotent, keep VCF headers clean and indexed, create output folders before writing, pin container versions, and keep LF line endings. Small habits, big difference.

