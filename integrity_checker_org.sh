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
#   • BAM   : check that the header is readable 
#           : check that the file ends with a valid BAM EOF marker (28 bytes with specific byte values)
#           : check if aligned (contains primary mapped alignments) vs unaligned
#           : check sortedness by coordinate (reported in header)
#           : check is file is human with refgenDetector 
#   • CRAM  : inspect file structure with samtools quickcheck
#           : check that the header is readable and contains @SQ lines  
#           : MD5 (M5) tags
#           : check if file is human with refgenDetector
#   • VCF   : verify that all samples in the VCF are present in the provided metadata CSV file           
#           : run VCFX_validator and capture errors
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

#function to initialize state for a new file; resets the _OKS, _FAILS, and _ERRS arrays to empty and sets the current file type and path for error reporting
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

# function to print the final status message for the current file based on the contents of the _OKS, _FAILS, and _ERRS arrays. 
# If there are any errors, it prints an [ERROR] message with all errors concatenated. If there are no errors but there are failures, it prints a [FAIL] message with all failures concatenated. If there are no errors or failures, it prints an [OK] message.
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

# help message function; prints usage instructions and exits with status 1
usage() {
    echo "Usage: $0 [OPTIONS] -f <file> [-s <samples>]"
    echo "Options:"
    echo "  -r, --run           Run mode (default)"
    echo "  -a, --analysis      Analysis mode (skips some checks)"
    echo "  -f, --file FILE     Input file to check"
    echo "  -s, --samples FILE    Metadata CSV file used for VCF sample-name checks"
    echo "  -h, --help          Show this help message and exit"
    exit 1
}

# helper function to control the mode errors
reject_file() {
    local type="$1" f="$2" msg="$3"
    begin_file "$type" "$f"
    err "$type" "$f" "$msg"
    end_file
    return $?
}

# helper function to control fails (e.g., missing helper tool)
internal_fail_file() {
    local type="$1" f="$2" msg="$3"
    begin_file "$type" "$f"
    fail "$type" "$f" "$msg"
    end_file
    return $?
}

FASTQ_LINES=40000      # FASTQ lines inspected
#BAM_EOF_BYTES=32768    # bytes read from end of BAM for EOF validation
#VCF_RECORDS=10000      # VCF/BCF records parsed
mode=""                # run or analysis
file=""                # input file path
samples=""             # metadata CSV file for VCF sample-name checks

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
    tmp=$(mktemp) || {
    fail "$type" "$f" "could not create temporary file"
    end_file; return $?
    }
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
        errs="fastQvalidator failed, but no detailed error message was returned. Please check if FASTQ format is valid"
    else
        errs="fastQvalidator failed in the first ${lines} lines: $errs"
    fi

    err "$type" "$f" "$errs"
    end_file; return $?
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
    mapped_primary=$(samtools view -h "$f" | head -n 100000 |samtools view -c -F 0x904 - 2>/dev/null)

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
        if ! command -v refgenDetector_main.py >/dev/null 2>&1; then
            fail "$type" "$f" "refgenDetector not found"
            end_file; return $?
        fi

        # Check if human 
        species=$(refgenDetector_main.py -f "$f" -t BAM/CRAM 2>/dev/null \
            | awk -F'Species detected:[[:space:]]*' '/Species detected:/ {print $2}' \
            | xargs)

        if [[ -z "$species" ]]; then
            err "$type" "$f" "refgenDetector produced no species result"
        elif [[ "$species" == "Homo sapiens" ]]; then
            ok "$type" "$f" "species: Homo sapiens"
        else
            err "$type" "$f" "refgenDetector: species is not human ($species)"
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

    local hdr_tmp sorted species qc_out qc_rc

    #samtools check 
    if ! command -v samtools >/dev/null 2>&1; then
        fail "$type" "$f" "samtools not found; CRAM check skipped"
        end_file; return $?
    fi

    hdr_tmp=$(mktemp) || {
        fail "$type" "$f" "could not create temporary file"
        end_file; return $?
    }
    trap 'rm -f "$hdr_tmp"' RETURN

    # header readable
    if ! samtools view -H "$f" > "$hdr_tmp" 2>/dev/null; then
        err "$type" "$f" "CRAM header is missing or unreadable. Please upload a valid CRAM file."
        end_file; return $?
    fi
    ok "$type" "$f" "header readable"

    qc_out=$(samtools quickcheck -vvv "$f" 2>&1)
    qc_rc=$?

    if (( qc_rc != 0 )); then
    err "$type" "$f" "CRAM file failed samtools quickcheck: $(printf '%s' "$qc_out" | tr '\n' '; ' | sed 's/; $//')"
    end_file; return $?
    fi

    # require @SQ
    if grep -q '^@SQ' "$hdr_tmp"; then
        ok "$type" "$f" "@SQ reference sequence entries present"
    else
        err "$type" "$f" "CRAM header is missing @SQ reference sequence entries. Please upload a valid CRAM file."
        end_file; return $?
    fi

    # Sortedness by coordinate (header-based; @HD may be missing)
    sorted=$(grep -m1 '^@HD' "$hdr_tmp" | grep -oE 'SO:[^[:space:]]+' | cut -d: -f2)

    if [[ "$sorted" == "coordinate" ]]; then
        ok "$type" "$f" "sortedness: coordinate"
    else
        err "$type" "$f" "CRAM not sorted by coordinate (SO:${sorted:-missing})"
    fi

    # require M5 on every @SQ line
    if awk 'BEGIN{ok=1} /^@SQ/ && $0 !~ /(^|[[:space:]])M5:/ {ok=0} END{exit ok?0:1}' "$hdr_tmp"; then
        ok "$type" "$f" "M5 reference MD5 tags present"
    else
        err "$type" "$f" "CRAM header is missing required reference MD5 (M5) tags. Please recreate the CRAM with the correct reference."
    fi

    # Human reference genome check
    if ! command -v refgenDetector_main.py >/dev/null 2>&1; then
        fail "$type" "$f" "refgenDetector not found"
        end_file; return $?
    fi

    species=$(
        refgenDetector_main.py -f "$f" -t BAM/CRAM 2>/dev/null \
        | awk -F'Species detected:[[:space:]]*' '/Species detected:/ {print $2}' \
        | xargs
    )

    if [[ -z "$species" ]]; then
        err "$type" "$f" "refgenDetector produced no species result"
    elif [[ "$species" == "Homo sapiens" ]]; then
        ok "$type" "$f" "species: Homo sapiens"
    else
        err "$type" "$f" "refgenDetector: species is not human ($species)"
    fi

    end_file; return $?
}

##############################################################################
# VCF
##############################################################################
check_vcf() {
    local f="$1" type="VCF"
    begin_file "$type" "$f"

    local vout rc errs missing
    local csv_tmp vcf_tmp

    # ----- VCFX_validator check -----
    if ! command -v VCFX_validator >/dev/null 2>&1; then
        fail "$type" "$f" "VCFX_validator not found"
        end_file; return $?
    fi

    # ----- bcftools check -----
    if ! command -v bcftools >/dev/null 2>&1; then
        fail "$type" "$f" "bcftools not found; sample name checks cannot be performed"
        end_file; return $?
    fi

    # ----- Check that all samples in the VCF are present in the metadata -----
    csv_tmp=$(mktemp) || {
        fail "$type" "$f" "Failed to create temporary file for CSV processing"
        end_file; return $?
    }
    trap 'rm -f "$csv_tmp" "$vcf_tmp"' RETURN

    vcf_tmp=$(mktemp) || {
        fail "$type" "$f" "Failed to create temporary file for VCF processing"
        end_file; return $?
    }

    # ----- input format / compression sanity check -----
    case "$f" in
        *.vcf.gz)
            if ! gzip -t "$f" 2>/dev/null; then
                err "$type" "$f" "Compressed VCF file could not be decompressed. Please check the file and resubmit."
                end_file; return $?
            fi
            ;;
        *.vcf.bz2)
            if ! bzip2 -t "$f" 2>/dev/null; then
                err "$type" "$f" "Compressed VCF file could not be decompressed. Please check the file and resubmit."
                end_file; return $?
            fi
            ;;
        *.vcf)
            ;;
    esac

    # Extract non-empty entries from first 3 CSV columns, skipping header
    awk -F',' 'NR > 1 {
        for (i = 1; i <= 3; i++) {
            gsub(/^[ \t]+|[ \t\r]+$/, "", $i)
            if ($i != "") print $i
        }
    }' "$samples" | sort -u > "$csv_tmp" || {
        fail "$type" "$f" "Could not read sample metadata file"
        end_file; return $?
    }

    # ----- extract sample names from VCF -----
    case "$f" in
        *.vcf.gz|*.vcf)
            if ! bcftools query -l "$f" > "$vcf_tmp" 2>/dev/null; then
                err "$type" "$f" "Could not read sample names from the VCF file. Please check the file and resubmit."
                end_file; return $?
            fi
            ;;
        *.vcf.bz2)
            if ! bzcat "$f" 2>/dev/null | bcftools query -l - > "$vcf_tmp" 2>/dev/null; then
                err "$type" "$f" "Could not read sample names from the VCF file. Please check the file and resubmit."
                end_file; return $?
            fi
            ;;
    esac

    missing=$(grep -Fxv -f "$csv_tmp" "$vcf_tmp" || true)
    if [[ -n "$missing" ]]; then
        err "$type" "$f" "Samples present in VCF but missing from registered metadata: $(printf '%s' "$missing" | tr '\n' ',' | sed 's/,$//')"
    fi

    # ----- Run VCFX_validator -----
    case "$f" in
        *.vcf.gz)
            vout=$(zcat "$f" 2>/dev/null | VCFX_validator 2>&1)
            rc=$?
            ;;
        *.vcf.bz2)
            vout=$(bzcat "$f" 2>/dev/null | VCFX_validator 2>&1)
            rc=$?
            ;;
        *.vcf)
            vout=$(VCFX_validator -i "$f" 2>&1)
            rc=$?
            ;;
    esac

    if printf '%s\n' "$vout" | grep -q '^Status:[[:space:]]*PASSED'; then
        ok "$type" "$f" "VCF file passed validation checks"
        end_file; return $?
    fi

    errs=$(
        printf '%s\n' "$vout" \
        | grep -E '^Error:' \
        | tr '\n' '; ' \
        | sed 's/; $//'
    )

    if [[ -z "$errs" ]]; then
        errs="VCF validation failed, but VCFX_validator did not return a detailed error message. Please check the file and consult the VCFX_validator documentation."
    else
        errs="VCFX_validator failed: $errs"
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
    fail "[FAIL] FILE $file - not found"
    exit 1
fi

# -- determine file type and execute checks -- 
case "$file" in
  *.fastq|*.fastq.gz|*.fq|*.fq.gz)
    if [[ "$mode" == "run" ]]; then
        check_fastq "$file"
        exit $?
    else
        reject_file "FASTQ" "$file" "FASTQ files need to be uploaded as RUNs"
        exit $?
    fi
    ;;

  *.bam|*.bam.gz)
    check_bam "$file"
    exit $?
    ;;
  *.cram|*.cram.gz) 
    if [[ "$mode" == "analysis" ]]; then
        check_cram "$file"
        exit $?
    else 
        reject_file "CRAM" "$file" "CRAM files need to be uploaded as ANALYSIS"
        exit $?
    fi
    ;;
    
  *.vcf|*.vcf.gz|*.vcf.bz2)
    if [[ "$mode" == "run" ]]; then
        reject_file "VCF" "$file" "VCF/BCF files need to be uploaded as ANALYSIS"
        exit $?
    fi
    
    if [[ -z "$samples" ]]; then
        internal_fail_file "VCF" "$file" "Sample metadata CSV (-s) was not provided to the QC script"
        exit $?
    fi

    check_vcf "$file"
    exit $?
    ;;
  *)
    echo "[WARNING] FILE $file - unsupported extension; skipping"
    exit 0
    ;;
esac

