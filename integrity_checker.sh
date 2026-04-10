#!/usr/bin/env bash
#
# 
# Purpose:
#   Perform lightweight integrity checks on common NGS file types before they
#   enter downstream pipelines. Exactly ONE line per file is printed:
#     [OK] / [FAIL] / [ERROR] <TYPE> <PATH> - <message>
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
#   • Captures all ERRORS and returns them.
#   • FAIL is non-fatal (e.g., missing helper tool); processing continues. Will be captured by dev team and logged in db
#   • OK indicates the file passed the implemented checks.
# ---------------------------------------------------------------------------

set +e  # aggregates errors; we do not want to abort afer the first failure but collect all possible issues
set -o pipefail # ensures that if any command in a pipeline fails, the entire pipeline is considered to have failed (i.e., it will return a non-zero exit status). This is important for error handling, especially when using tools like grep or awk in pipelines, as it allows us to detect failures that might otherwise be masked by successful commands later in the pipeline.

declare -a _OKS _FAILS _ERRS
_cur_type=""
_cur_file=""

begin_file() {
    _cur_type="$1"
    _cur_file="$2"
    _OKS=()
    _FAILS=()
    _ERRS=()
}

ok()   { _OKS+=("$3"); }
fail() { _FAILS+=("$3"); }
err()  { _ERRS+=("$3"); }


end_file() {
  local msg

  if ((${#_ERRS[@]})); then
    msg=$(printf '%s; ' "${_ERRS[@]}")
    msg=${msg%; }   # remove trailing "; "
    printf '[ERROR] %s %s - %s\n' "$_cur_type" "$_cur_file" "$msg"
    return 1
  fi

  if ((${#_FAILS[@]})); then
    msg=$(printf '%s; ' "${_FAILS[@]}")
    msg=${msg%; }
    printf '[FAIL] %s %s - %s\n' "$_cur_type" "$_cur_file" "$msg"
    return 0
  fi

  printf '[OK] %s %s - file OK\n' "$_cur_type" "$_cur_file"
  return 0
}

usage() {
    echo "Usage: $0 [OPTIONS] -f <file> -s <samples>"
    echo "Options:"
    echo "  -r, --run           Run mode (default)"
    echo "  -a, --analysis      Analysis mode (skips some checks)"
    echo "  -f, --file FILE     Input file to check"
    echo "  -s, --samples SAMPLES   Sample identifiers (comma-separated)"
    echo "  -h, --help          Show this help message and exit"
    exit 1
}

#FASTQ_LINES=40000      # FASTQ lines inspected
BAM_EOF_BYTES=32768    # bytes read from end of BAM for EOF validation
VCF_RECORDS=10000      # VCF/BCF records parsed
mode=""                # run or analysis
file=""                # input file path
samples=""             # sample identifiers (comma-separated)

#get options and arguments 
while [[ $# -gt 0 ]]; do
    case "$1" in 
        -r|--run) mode="run"; shift ;; #shift removes this option from the list of arguments; next iteration will process the next one
        -a|--analysis) mode="analysis"; shift ;;
        -f|--file) 
            [[ $# -lt 2 ]] && { echo "ERROR: -file requires an argument"; usage;}
            file="$2"; shift 2 ;; 
        -s|--samples) 
            [[ $# -lt 2 ]] && { echo "ERROR: -samples requires an argument"; usage;}
            samples="$2"; shift 2 ;;
        -h|--help) usage ;;
        *) echo "ERROR: Unknown option $1"; usage ;;
    esac
done

##############################################################################
# FASTQ
##############################################################################
check_fastq() {
    local f="$1" type="FASTQ" lines=$(( FASTQ_LINES - (FASTQ_LINES % 4) ))
    begin_file "$type" "$f"

    # check if file is empty
    [[ -s "$f" ]] || {
        err "$type" "$f" "File is empty. Please upload non-empty file"
        end_file; return $?
    }

    local tmp
    tmp=$(mktemp)
    trap 'rm -f "$tmp"' RETURN

    if [[ "$f" == *.gz ]]; then
        # verify gzip integrity (reads fill compressed stream, no output written)
        if ! gzip -t "$f" 2>/dev/null; then
            err "$type" "$f" "Compressed file is corrupted/incomplete. Re-create .gz file and upload again."
            end_file; return $?
        fi

        # extract only first sample chhunk; pipefall is disabeled here to prevent 
        # head exits early by design and can trigger SIGPIPE in gunzip, which would cause the entire pipeline to fail if pipefail were enabled. 
        if ! (
            set +o pipefail
            gunzip -c "$f" 2>/dev/null | head -n "$lines" > "$tmp"
        ); then  
            err "$type" "$f" "Failed to read compressed file. Please check file integrity and re-upload."
            end_file; return $?
        fi
    else
        # for uncompressed files, read only sample prefix
        if ! head -n "$lines" "$f" > "$tmp" 2>/dev/null; then
            err "$type" "$f" "Failed to read file. Please check file integrity and re-upload."
            end_file; return $?
        fi
    fi

    # basic format check: first line should start with '@'
    if ! head -n 1 "$tmp" | grep -q '^@'; then
    err "$type" "$f" "Format error: record 1 does not start with '@'. Please upload valid FASTQ file."
    end_file; return $?
    fi

    # validator check
    if ! command -v fastQValidator >/dev/null 2>&1; then
        fail "$type" "$f" "fastQValidator not found"
        end_file; return $?
    fi

    local vout rc errs
    vout=$(fastQValidator --file "$tmp" --disableSeqIDCheck 2>&1) #runs fastQValidator on extracted sample chunk
    rc=$?

    if (( rc == 0 )); then
        ok "$type" "$f" "validator passed"
        end_file; return $?
    fi

    errs=$(
        printf '%s\n' "$vout" \
        | grep '^ERROR on Line ' \
        | tr '\n' '; ' \
        | sed 's/; $//'
    )

    if [[ -z "$errs" ]]; then 
        errs="FASTQ format validation failed, but no detailed error message was returned. Please check if FASTQ format is valid"
    else
        errs="FASTQ format validation failed in the first ${lines} lines: $errs"
    fi

    err "$type" "$f" "$errs"
    end_file; return $?
    fi
}

##############################################################################
# BAM
##############################################################################
check_bam() {
    local f="$1" type="BAM"
    begin_file "$type" "$f"

    local hdr_tmp sorted species align="False" mapped_primary
    local eof_hex expected_eof="1f8b08040000000000ff0600424302001b0003000000000000000000"
    
    # ----- samtools check -----
    if ! command -v samtools >/dev/null 2>&1; then
        fail "$type" "$f" "samtools not found; header/sortedness checks skipped"
        end_file
        return $?
    fi

    # initilize empty hdr_temp and trap for cleanup
    hdr_tmp=$(mktemp) || {
    fail "$type" "$f" "could not create temporary file"
    end_file; return $?
    }
    trap 'rm -f "$hdr_tmp"' RETURN

    #check if header is readable
    if ! samtools view -H "$f" > "$hdr_tmp" 2>/dev/null; then
        err "$type" "$f" "BAM header missing or unreadable"
        end_file
        return $?      
    fi
    ok "$type" "$f" "header readable"

    # Check EOF marker
    eof_hex=$(tail -c 28 "$f" 2>/dev/null | xxd -p -c 28)

    if [[ "$eof_hex" == "$expected_eof" ]]; then
        ok "$type" "$f" "BAM EOF marker present"
    else
        err "$type" "$f" "BAM EOF marker missing or file truncated. Please recreate the BAM file and resubmit."
    fi

    # Check if BAM is aligned or unaligned 
    mapped_primary=$(samtools view -c -F 0x904 "$f" 2>/dev/null)

    #check if mapped_primary is empty (due to error)
    if [[ -z "$mapped_primary" ]]; then
    err "$type" "$f" "Could not inspect BAM alignment records. Please check the file and resubmit."
    end_file; return $?
    fi

    # alignment check
    if (( mapped_primary > 0)); then
        align="True"
        if [[ "$mode" == "run" ]]; then
            err "$type" "$f" "This BAM appears to be ALIGNED (contains primary mapped alignments); Please upload this file as an ANALYSIS, not a RUN"
        else 
            ok "$type" "$f" "BAM appears to be ALIGNED + uploaded as ANALYSIS"
        fi
    else
        align="False"
        if [[ "$mode" == "run" ]]; then
            ok "$type" "$f" "BAM appears to be UNALIGNED (does not contain primary mapped alignments) + uploaded as RUN"
        else 
            err "$type" "$f" "This BAM appears to be UNALIGNED (does not contain primary mapped alignments); Please upload this file as a RUN, not an ANALYSIS"
        fi
    fi

    # for aligned files 
    if [[ "$align" == "True" ]]; then
        # check sortedness by coordinate (reported in header)
        sorted=$(grep -m1 '^@HD' "$hdr_tmp" | grep -oE "SO:[^[:space:]]*" | cut -d: -f2)
        if [[ "$sorted" == "coordinate" ]] ; then
            ok "$type" "$f" "BAM file sorted by coordinate"
        else
            err "$type" "$f" "BAM not sorted by coordinate (SO:${sorted:-missing})"
        fi

        # Check if refgenDetector is available
        if ! command -v refgenDetector >/dev/null 2>&1; then
            fail "$type" "$f" "refgenDetector not found"
            end_file; return $?
        fi

        # Check if human 
        species=$(refgenDetector -f "$f" -t BAM/CRAM 2>/dev/null \
            | awk -F'Species detected:[[:space:]]*' '/Species detected:/ {print $2}' \
            | xargs)

        if [[ -z "$species" ]]; then
            err "$type" "$f" "refgenDetector produced no species result"
        elif [[ "$species" == "Homo sapiens" ]]; then
            ok "$type" "$f" "species: Homo sapiens"
        else
            err "$type" "$f" "species is not human ($species)"
        fi
    fi
        end_file; return $?

}

##############################################################################
# CRAM
##############################################################################
check_cram() {
    local f="$1" type="CRAM"
    begin_file "$type" "$f"

    local hdr_tmp sorted species

    #samtools check 
    if ! command -v samtools >/dev/null 2>&1; then
        fail "$type" "$f" "samtools not found; CRAM check skipped"
        end_file; return $?
    fi

    # Read header once 
    hdr_tmp=$(mktemp)
    if ! samtools view -H "$f" > "$hdr_tmp" 2>/dev/null; then
        err "$type" "$f" "CRAM header missing or unreadable"
        rm -f "$hdr_tmp"
        end_file; return $?
    fi

    # Require @SQ (reference dictionary)
    if grep -q '^@SQ' "$hdr_tmp"; then
        ok "$type" "$f" "header readable; @SQ present"
    else
        err  "$type" "$f" "header missing @SQ"
        rm -f "$hdr_tmp"
        end_file; return $?
    fi

    # Sortedness by coordinate (header-based; @HD may be missing)
    sorted=$(grep -m1 '^@HD' "$hdr_tmp" | grep -oE 'SO:[^[:space:]]+' | cut -d: -f2)

    if [[ "$sorted" == "coordinate" ]]; then
        ok "$type" "$f" "sortedness: coordinate"
    else
        err "$type" "$f" "CRAM not sorted by coordinate (SO:${sorted:-missing})"
    fi

    # Reference MD5 tags (now REQUIRED)
    if grep -q 'M5:' "$hdr_tmp"; then
        ok "$type" "$f" "M5 reference MD5 tags present"
    else
        err "$type" "$f" "missing required M5 reference MD5 tags"
    fi

    rm -f "$hdr_tmp"

    # Human reference genome check
    if ! command -v refgenDetector >/dev/null 2>&1; then
        fail "$type" "$f" "refgenDetector not found"
        end_file; return $?
    fi

    species=$(
        refgenDetector -f "$f" -t BAM/CRAM 2>/dev/null \
        | awk -F'Species detected:[[:space:]]*' '/Species detected:/ {print $2}' \
        | xargs
    )

    if [[ -z "$species" ]]; then
        err "$type" "$f" "refgenDetector produced no species result"
    elif [[ "$species" == "Homo sapiens" ]]; then
        ok "$type" "$f" "species: Homo sapiens"
    else
        err "$type" "$f" "species is not human ($species)"
    fi

    end_file; return $?
}

##############################################################################
# VCF / BCF
##############################################################################
check_vcf() {
    local f="$1" type="VCF/BCF"
    begin_file "$type" "$f"

    # ----- tool check -----
    if ! command -v VCFX_validator >/dev/null 2>&1; then
        fail "$type" "$f" "VCFX_validator not found"
        end_file; return $?
    fi

    local vout rc errs

    if [[ "$f" == *.gz ]]; then
        vout=$(zcat "$f" | VCFX_validator 2>&1)
        rc=$?
    elif [[ "$f" == *.bz2 ]]; then
        vout=$(bzcat "$f" | VCFX_validator 2>&1)
        rc=$?
    else
        vout=$(VCFX_validator -i "$f" 2>&1)
        rc=$?
    fi

    # PASS condition: explicit status line
    if printf '%s\n' "$vout" | grep -q '^Status:[[:space:]]*PASSED'; then
        ok "$type" "$f" "validator passed"
        end_file; return $?
    fi

    # Collect all error lines (join with '; ' like FASTQ)
    errs=$(
        printf '%s\n' "$vout" \
        | grep -E '^Error:' \
        | tr '\n' '; ' \
        | sed 's/; $//'
    )

    if [[ -z "$errs" ]]; then
        # fallback: if validator failed but didn't print "Error:" lines
        errs=$(printf '%s\n' "$vout" | head -n 1)
        [[ -z "$errs" ]] && errs="VCFX_validator failed (rc=$rc) with no output"
    fi

    err "$type" "$f" "$errs"
    end_file; return $?
}


##############################################################################
# Main
##############################################################################

# --validate required arguments --
if [[ -z "$file" ]]; then
  echo "Error: -file is required"
  usage
fi

if [[ -z "$mode" ]]; then
  echo "Error: you must specify either -run(-r) or -analysis(-a)"
  usage
fi

# check file exists
if [[ ! -f "$file" ]]; then
    echo "[ERROR] FILE $file - not found"
    exit 1
fi

# if [[ -z "$samples" ]]; then
#   echo "Error: -samples is required"
#   usage
# fi

# -- determine file type and execute checks -- 
case "$file" in
  *.fastq|*.fastq.gz|*.fq|*.fq.gz)
    if [[ "$mode" == "run" ]]; then
        check_fastq "$file"
        exit $?
    fi
    if [[ "$mode" == "analysis" ]]; then
        err "FASTQ" "$file" "FASTQ files need to be uploaded as RUNs" 
    fi
    ;;
  *.bam|*.bam.gz)
    if [[ "$mode" == "run" ]]; then
        check_unaligned_bam "$file"
        exit $?
    fi
    if [[ "$mode" == "analysis" ]]; then
        check_aligned_bam "$file"
        exit $?
    fi
    ;;
  *.cram|*.cram.gz) # @@@ TODO: discuss if CRAM files also can be alsigned/unaligned and define checks to be performed
    check_cram "$file"
    exit $?
    ;;
  *.vcf|*.vcf.gz|*.bcf|*.bcf.gz|*.vcf.bz2|*.bcf.bz2)
    
    if [[ "$mode" == "run" ]]; then
        err "VCF/BCF" "$file" "VCF/BCF files need to be uploaded as ANALYSIS"
    fi
    if [[ "$mode" == "analysis" ]]; then
        check_vcf "$file"
        exit $?
    fi
    ;;
  *)
    echo "[WARNING] FILE $file - unsupported extension; skipping"
    exit 0
    ;;
esac

