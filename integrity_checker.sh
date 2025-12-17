#!/usr/bin/env bash
#
# 
# Purpose:
#   Perform lightweight integrity checks on common NGS file types before they
#   enter downstream pipelines. Exactly ONE line per file is printed:
#     [OK] / [WARNING] / [ERROR] <TYPE> <PATH> - <message>
#
# Checks implemented
#   • FASTQ : extract the first 40 000 lines and run fastQValidator
#   • BAM   : inspect the first 500 header lines and verify the BAM EOF marker
#   • CRAM  : inspect the first 500 header lines and REQUIRE reference
#             MD5 (M5) tags (missing M5 ⇒ ERROR)
#   • VCF   : parse the header plus the first 10 000 variant records with
#             bcftools head (fatal parse ⇒ ERROR)
#           : verify that VCF/BCF files are sorted according to HTSlib rules
#
# Exit status
#   • Any ERROR terminates execution immediately (set -e).
#   • WARNING is non-fatal (e.g., missing helper tool); processing continues.
#   • OK indicates the file passed the implemented checks.
# ---------------------------------------------------------------------------

set -e  # abort on command failure

FASTQ_LINES=40000      # FASTQ lines inspected
BAM_EOF_BYTES=32768    # bytes read from end of BAM for EOF validation
VCF_RECORDS=10000      # VCF/BCF records parsed
HUM_REF_GENOME=("hg16" "hg17" "hg18" "GRCh37" "GRCh38" "T2T")  # human reference genome names

ok()   { echo "[OK] $1 $2 - $3"; }
fail() { echo "[FAIL] $1 $2 - $3"; }
err()  { echo "[ERROR] $1 $2 - $3"; exit 1; }


##############################################################################
# FASTQ
##############################################################################
check_fastq() {
    local f="$1" type="FASTQ"
    local tmp; tmp=$(mktemp)

    case "$f" in
        *.gz)  gunzip -c "$f" | head -n "$FASTQ_LINES" > "$tmp" ;;
        *)     head    -n "$FASTQ_LINES" "$f"         > "$tmp" ;;
    esac

    if ! command -v fastQValidator >/dev/null 2>&1; then
        rm -f "$tmp"
        fail "$type" "$f" "fastQValidator not found; FASTQ check skipped"
        return
    fi

    if fastQValidator --file "$tmp" --maxErrors 1 --disableSeqIDCheck >/dev/null 2>&1; then
        rm -f "$tmp"
        ok "$type" "$f" "validator passed on first ${FASTQ_LINES} lines"
    else
        rm -f "$tmp"
        err "$type" "$f" "fastQValidator reported format problems"
    fi
}

##############################################################################
# BAM
##############################################################################
check_bam() {
    local f="$1" type="BAM"

    if ! command -v samtools >/dev/null 2>&1; then
        fail "$type" "$f" "samtools not found; BAM check skipped"
        return
    fi

    # Header presence
    if ! samtools view -H "$f" 2>/dev/null | head -n 500 | grep -Eq '^@HD|^@SQ'; then
        err "$type" "$f" "header missing @HD/@SQ"
    fi

    # Sortedness by coordinate
    sorted=$(samtools view -H "$f" | grep "@HD" | grep -o "SO:[^[:space:]]*" | cut -d: -f2)
    if [[ "$sorted" == "coordinate" ]] ; then
        ok "$type" "$f" "header OK; @CO lines present"
    else
        err "$type" "$f" "BAM not sorted by coordinate"
    fi

    # EOF marker
    local tail_tmp; tail_tmp=$(mktemp)
    tail -c "$BAM_EOF_BYTES" "$f" > "$tail_tmp" 2>/dev/null || true
    if ! xxd -p "$tail_tmp" | grep -iq '42430200'; then
        rm -f "$tail_tmp"
        err "$type" "$f" "BAM EOF magic (42430200) not found"
    fi
    rm -f "$tail_tmp"

    ok "$type" "$f" "header OK; EOF magic present"

    # Human reference genome check
    if ! command -v refgenDetector_main.py >/dev/null 2>&1; then
        fail "$type" "$f" "refgenDetector not found; human reference genome check skipped"
        return
    fi

    ref=$(refgenDetector_main.py -f "$f" -t BAM/CRAM 
        | awk -F'Species detected:[[:space:]]*' '/Species detected:/ {print $2}' \
        | xargs)

    if [[ "$ref" == "Homo sapiens" ]]; then
        ok "$type" "$f" "BAM file is human"
    else
        err "$type" "$f" "BAM file is not human"
    fi

}


##############################################################################
# CRAM
##############################################################################
check_cram() {
    local f="$1" type="CRAM"

    if ! command -v samtools >/dev/null 2>&1; then
        fail "$type" "$f" "samtools not found; CRAM check skipped"
        return
    fi

    # Header presence
    if ! samtools view -H "$f" 2>/dev/null | head -n 500 | grep -Eq '^@HD|^@SQ'; then
        err "$type" "$f" "header missing @HD/@SQ"
    fi

    # Sortedness by coordinate
    sorted=$(samtools view -H "$f" | grep "@HD" | grep -o "SO:[^[:space:]]*" | cut -d: -f2)
    if [[ "$sorted" == "coordinate" ]] ; then
        ok "$type" "$f" "header OK; @CO lines present"
    else
        err "$type" "$f" "CRAM not sorted by coordinate"
    fi

    # Reference MD5 tags (now REQUIRED)
    if samtools view -H "$f" 2>/dev/null | grep -q 'M5:'; then
        ok "$type" "$f" "header OK; M5 reference MD5 tags present"
    else
        err "$type" "$f" "header OK; missing required M5 reference MD5 tags"
    fi

    # Human reference genome check
    if ! command -v refgenDetector_main.py >/dev/null 2>&1; then
        fail "$type" "$f" "refgenDetector not found; human reference genome check skipped"
        return
    fi

    ref=$(refgenDetector_main.py -f "$f" -t BAM/CRAM 
        | awk -F'Species detected:[[:space:]]*' '/Species detected:/ {print $2}' \
        | xargs)

    if [[ "$ref" == "Homo sapiens" ]]; then
        ok "$type" "$f" "CRAM file is human"
    else
        err "$type" "$f" "CRAM file is not human"
    fi
}

##############################################################################
# VCF / BCF
##############################################################################
check_vcf() {
    local f="$1" type="VCF/BCF"

    if ! command -v bcftools >/dev/null 2>&1; then
        fail "$type" "$f" "bcftools not found; VCF/BCF check skipped"
        return
    fi

    if bcftools head -n "$VCF_RECORDS" "$f" >/dev/null 2>&1; then
        ok "$type" "$f" "bcftools parsed header + first ${VCF_RECORDS} records"
    else
        err "$type" "$f" "bcftools failed parsing within first ${VCF_RECORDS} records"
    fi
}

check_variant_sorted() {
    local f="$1" type="VCF/BCF"

    if [[ -z "$f" ]]; then
        fail "$type" "$f" "no file given; sortedness check skipped"
        return
    fi

    # If BCF-like, we need bcftools to convert to VCF text
    if [[ "$f" == *.bcf || "$f" == *.bcf.gz || "$f" == *.bcf.bz2 ]]; then
        if ! command -v bcftools >/dev/null 2>&1; then
            fail "$type" "$f" "bcftools not found; sortedness check skipped for BCF"
            return
        fi
    fi

    # Run the sortedness check; capture first offending locus if unsorted.
    # awk prints "chr:pos" and exits 1 on the first out-of-order record,
    # or prints nothing and exits 0 if everything is sorted.
    local unsorted_at=""
    if ! unsorted_at=$(
        if [[ "$f" == *.bcf || "$f" == *.bcf.gz || "$f" == *.bcf.bz2 ]]; then
            # Any BCF (compressed or not) → VCF stream via bcftools
            bcftools view -Ov "$f" 2>/dev/null
        elif [[ "$f" == *.vcf.bz2 ]]; then
            bzcat "$f" 2>/dev/null
        elif [[ "$f" == *.vcf.gz ]]; then
            zcat "$f" 2>/dev/null
        else
            cat "$f"
        fi | awk -F'\t' '
            BEGIN {
                last_tid = -1
                last_pos = -1
                nctg = 0
            }

            # Collect contigs from header
            /^##contig=/ {
                if (match($0, /ID=([^,>]+)/, a)) {
                    cid = a[1]
                    if (!(cid in tid)) {
                        tid[cid] = nctg
                        nctg++
                    }
                }
                next
            }

            # Skip other header lines
            /^#/ { next }

            # Variant lines
            {
                chrom = $1
                pos   = $2 + 0

                # If chrom not in header, append at end
                if (!(chrom in tid)) {
                    tid[chrom] = nctg
                    nctg++
                }
                t = tid[chrom]

                # Same logic as HTSlib: (tid, pos) must not go backwards
                if (t < last_tid || (t == last_tid && pos < last_pos)) {
                    # Print first offending locus to stdout and exit non-zero
                    printf("%s:%d\n", chrom, pos)
                    exit 1
                }

                last_tid = t
                last_pos = pos
            }
        '
    ); then
        # Non-zero exit from awk / pipeline
        if [[ -n "$unsorted_at" ]]; then
            err "$type" "$f" "variants not sorted by contig/POS (first offending locus: $unsorted_at)"
        else
            err "$type" "$f" "sortedness scan failed (parse or I/O error)"
        fi
    else
        ok "$type" "$f" "variants sorted by contig/POS"
    fi
}





##############################################################################
# Main
##############################################################################
for file in "$@"; do
    if [[ ! -f "$file" ]]; then
        fail "FILE" "$file" "not found; skipping"
        continue
    fi

    case "$file" in
        *.fastq|*.fastq.gz|*.fq|*.fq.gz)
            check_fastq "$file" ;;
        *.bam|*.bam.gz)
            check_bam "$file" ;;
        *.cram|*.cram.gz)
            check_cram "$file" ;;
        *.vcf|*.vcf.gz|*.bcf|*.bcf.gz|*.vcf.bz2|*.bcf.bz2)
            check_vcf "$file" 
            check_variant_sorted "$file" ;; 
        *)
            fail "FILE" "$file" "unsupported extension; skipping" ;;
    esac
done

