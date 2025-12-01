# Portable Precision – GIAB Parity (Nextflow DSL2)

Hello, This is  **Peter Oluwatimilehin** 👋  

This repo is my flagship personal project as I move from:

- a **Microbiology** background,
- hands-on work running an **electronic health record (EHR)** at a busy student clinic,
- and an online **cancer genomics scholarship**,

into a career in **health data science / bioinformatics**.

I wanted something practical that lives in the same world as modern health data work:  
small enough to run on a laptop, but serious enough to look like a real germline pipeline that might one day sit close to a trusted research environment or diagnostic workflow.

This is that project.

---

## 1. What this project is

**Portable Precision – GIAB Parity** is an end-to-end **germline variant-calling pipeline** written in **Nextflow DSL2**.

It’s built around the same ideas used in Genome in a Bottle (GIAB) benchmarking:

- start with FASTQ reads and a reference genome,
- run through alignment and variant calling,
- compare your calls to a **truth set** inside a defined “confident region” to see how well you’re doing.

Today, the pipeline focuses on:

- **Reproducibility** – everything is expressed as a Nextflow workflow, not a loose list of commands.
- **Portability** – the same logic can run:
  - locally on my machine,
  - with `workDir` and `outdir` on **S3**, and
  - in **AWS Batch** environments (EC2 and Fargate) using containers and IAM.

It’s not a production pipeline; it’s a **learning + portfolio pipeline** that is honest about what works, what’s still fragile, and what I’m learning in the process.

---

## 2. Why I built this (and why now)

This project wasn’t requested by any course, employer or programme.  
I built it because of how my experiences have layered over time:

- During my **BSc in Microbiology**, I learned how infections, immunity and mutations play out in real people, long before any data hits a spreadsheet.
- At **UHS Jaja Clinic** (University of Ibadan), I worked in the records office and eventually led the digitisation of more than 2,500 student records into an EHR, training colleagues on using SmartClinic and seeing directly how data quality affects care.
- In 2024, I joined an online **cancer genomics research scholarship**, focusing on how **non-synonymous SNPs in protein-coding genes** may influence cancer susceptibility and progression. That’s where I started reading about GIAB, variant benchmarking and tools like **hap.py** that compare callsets against curated truth datasets.
- Now, on my **MSc in Data Science** in the UK, I’m formalising the maths and methods behind what I’ve been feeling intuitively: statistical thinking, predictive modelling, and careful interpretation of metrics.

At some point it clicked:

> I don’t just want to run “black box” tools.  
> I want to understand what happens between raw FASTQs and a table of variants that clinicians and researchers can trust.

So I gave myself a project that:

1. **Feels like real work**  
   A pipeline that could plausibly be extended into something used in research or clinical benchmarking, not just a toy script.

2. **Forces me to think about portability**  
   Not only “does it run on my laptop?”, but “how does this behave with object storage, Batch schedulers, containers and IAM in the loop?”.

3. **Produces evidence, not just claims**  
   MultiQC reports, hap.py-style summary metrics, and Nextflow reports/timelines that you can actually open and inspect.

HDR UK and similar organisations are examples of places where this mindset matters, but this project exists because it aligns with where I want my career to go **beyond any single internship**.

---

## 3. What the pipeline currently does

Right now, the pipeline:

- Takes **paired-end FASTQ** files.
- Uses a small reference genome and (optionally) a GIAB-style **truth VCF + confident regions BED**.
- Runs through:

  1. **FastQC** – quality control on raw reads.  
  2. **BWA-MEM** – alignment to the reference genome.  
  3. **samtools** – sort, index, and generate basic alignment stats.  
  4. **GATK HaplotypeCaller** – germline variant calling.  
  5. **bcftools** – basic VCF statistics and filtering.  
  6. **MultiQC** – one HTML report that aggregates QC and stats for quick review.

Earlier versions also integrated **hap.py** to compare the called VCF against a truth set within a confident region, producing precision/recall/F1 and error summaries. Example outputs from that stage are still kept separately as evidence and are part of the roadmap to restore as a clean, first-class step.

Under the hood, it’s written in **Nextflow DSL2**, with modular processes and profiles for local and AWS Batch execution.

---

## 4. How it’s built (Nextflow DSL2)

The pipeline is written in **Nextflow DSL2** with one module per tool.

Key design points:

- **Symmetric A/B lanes** – mirroring GIAB-style comparisons (e.g. two callsets against the same truth).
- **Modular processes** – each module (`FastQC`, `BWA`, `samtools`, `GATK HC`, `bcftools`, `MultiQC`) can be swapped or reused.
- **Config-driven profiles** – a `local` profile for laptop runs and an `awsbatch_fg` profile for AWS Batch Fargate with S3-backed `workDir` and `outdir`.

Some of the things I’ve learned on the way:

- Channels in DSL2 are **single-consumer** – if two steps need the same data, you must split it correctly.
- If a module is used twice, **alias it** (e.g. `GATK_HC_A` and `GATK_HC_B`).
- Make reference/index steps **idempotent** so there are no file name collisions.
- `hap.py` needs clean, indexed VCFs with matching contigs/headers and a proper output folder.
- Pinning container versions and using LF line endings helps avoid “works on my machine” problems.

---

## 5. Cloud & portability (AWS – honest summary)

I’ve used this project to explore how a genomics workflow behaves off my laptop:

- Configured **S3** buckets so the pipeline can:
  - run locally with `workDir` and `outdir` on S3, and
  - mimic what it might feel like inside a TRE with object storage.
- Created **AWS Batch** compute environments and queues for:
  - EC2-backed jobs, and
  - Fargate-backed “serverless” jobs.
- Successfully ran **smaller Nextflow workflows** end-to-end on Batch.
- Ran the **core of this pipeline** (FastQC → BWA → samtools → GATK HC → bcftools) on Fargate with S3-backed storage on a separate branch.
- Debugged failures around MultiQC and hap.py on Fargate (container size, tools, permissions), and used that to sharpen my understanding of IAM roles, container images and Batch logs.

Some stages (especially MultiQC and hap.py on Fargate) are still work-in-progress.  
I’m keeping that visible on purpose – this repo is a record of **real learning**, not a fake “everything is perfect” story.

---

## 6. Quick start – tiny local run

This is intentionally tiny so it doesn’t overwhelm your machine.

### Requirements

- Linux / WSL2 or macOS  
- Conda / mamba  
- Nextflow (developed with v25.x)

### Setup

```bash
# clone
git clone https://github.com/PeterHDS/peterolu-giab-parity-nextflow.git
cd peterolu-giab-parity-nextflow

# create environment
conda create -n giab-nf -y \
  -c conda-forge -c bioconda \
  nextflow \
  bwa \
  samtools \
  fastqc \
  gatk4 \
  bcftools \
  multiqc

conda activate giab-nf

