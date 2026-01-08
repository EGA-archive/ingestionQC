# Multi‑format File Integrity Checker

This Bash script provides lightweight and rapid integrity checks for common next-generation sequencing (NGS) data files, enabling early detection of problematic or corrupted files prior to downstream bioinformatics analysis.

## Purpose

- Quickly identify corrupt or improperly formatted genomic data files.
- Ensure basic file format compliance, reducing the risk of pipeline failure.

## How it works

1. Detects file type from the extension.
2. Runs a lightweight check tailored to that format:

   * **FASTQ** → inspect first 40 000 lines → `fastQValidator`
   * **BAM** → header sanity check + EOF signature
   * **CRAM** → header sanity check + MD5 presence
   * **VCF/BCF** → header + first 10 000 records via `bcftools`
  
Each file produces a single summary line of output.

## Supported file types

### 1. FASTQ (`.fastq`, `.fastq.gz`, `.fq`, `.fq.gz`)

- Extracts the first **40,000 lines**.
- Performs validation using `fastQValidator` to detect critical formatting errors.
- Output:
  - `[OK]` if validator passes
  - `[ERROR]` if validator detects problems
  - `[WARNING]` if `fastQValidator` is missing (check skipped)

### 2. BAM (`.bam`, `.bam.gz`)

- Inspects the first **500 header lines** using `samtools` to confirm presence of essential tags (`@HD` and/or `@SQ`).
- Verifies presence of the BAM-specific EOF marker (`42430200`) within the last 32 KB, ensuring the file isn't truncated.
- Checks if file is sorted by corrdinates
- Checks if file is human with refgenDetector 
- Output:
  - `[OK]` if all checks are okay
  - `[ERROR]` collected error messages of all checks performed
  - `[WARNING]` if any tools are missing

### 3. CRAM (`.cram`, `.cram.gz`)

- Checks the first **500 header lines** using `samtools` for essential header tags (`@HD` and/or `@SQ`).
- Verifies presence of reference **MD5 checksum (M5)** tags in the header.
- Checks if file is sorted by corrdinates
- Checks if file is human with refgenDetector 
- Output:
  - `[OK]` if all checks are okay
  - `[ERROR]` collected error messages of all checks performed
  - `[WARNING]` if any tools are missing

### 4. VCF/BCF (`.vcf`, `.vcf.gz`, `.bcf`, `.bcf.gz`, `.vcf.bz2`, `.bcf.bz2`)

- Verifies compliance to VCF 4.x specifications and sortendess by VCFX_validator
- Rapidly identifies fatal syntax, format, or specification errors.
- Output:
  - `[OK]` if VCFX_validator runs without errors
  - `[ERROR]` VCFX_validator error message
  - `[WARNING]` if tool is missing

## Dependencies

| Tool                                                          | Purpose                |
| ------------------------------------------------------------- | ---------------------- |
| **bash** ≥4                                                   | scripting language     |
| GNU **coreutils** (`mktemp`,`head`, `tail`, `grep`, `xxd`, `awk`, `sed`, `cut`, `tr`, `xargs`, `gunzip`, `zcat`, `bzcat`)             | basic ops              |
| [`fastQValidator`](https://github.com/statgen/fastqvalidator) | FASTQ validation       |
| [`samtools`](https://www.htslib.org/)                         | BAM/CRAM header checks |
| [`refgenDetector`](https://github.com/EGA-archive/refgenDetector) | human check|
| [`VCFX_validator`](https://ieeta-pt.github.io/VCFX/installation/) | VCF/BCF validation |


Ensure each tool is on your `$PATH`.

## Usage

```bash
chmod +x integrity_checker.sh
./integrity_checker.sh file1.fastq.gz file2.bam file3.vcf.gz
```

The script auto‑detects the format and prints status messages for each file.

### Sample output

```text
[OK] FASTQ file1.fastq.gz - validator passed on first 40000 lines
[OK] BAM file2.bam - header OK; EOF magic present
[OK] VCF/BCF file3.vcf.gz - bcftools parsed header + first 10000 records
```

## Exit status and log messaging

| Status       | Description                                                             |
| ------------ | ----------------------------------------------------------------------- |
| `[OK]`       | File passed all checks.                                                 |
| `[WARNING]`  | Recommended software not available.                                     |
| `[ERROR]`    | Fatal issue; Errors are collected and printed separated by ; (one line).|

A pass indicates that the *sampled portion* is valid; it does **not** guarantee that the entire file is error‑free.

##  Limitations

* Only partial validation for speed, deep‑file corruption may go unnoticed.
* Default thresholds can be adjusted in the script:
* Validates only selected sections of the file:
  ```text
   *40k FASTQ lines
   *500 BAM/CRAM header lines
   *10k VCF records
```text
