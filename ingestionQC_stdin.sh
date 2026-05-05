#!/usr/bin/env bash

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

ok()   { _OKS+=("$1"); }
fail() { _FAILS+=("$1"); }
err()  { _ERRS+=("$1"); }

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
    echo "Usage: $0 [OPTIONS] -f <stdin> -e <extension> [-s <samples>]"
    echo "Options:"
    echo "  -e, --extension EXT      File extension (e.g., fastq.gz, bam, cram)"
    echo "  -r, --run           Run mode (default)"
    echo "  -a, --analysis      Analysis mode (skips some checks)"
    echo "  -f, --file STDIN     Input STDIN to check"
    echo "  -s, --samples FILE    Metadata CSV file used for VCF sample-name checks"
    echo "  -h, --help          Show this help message and exit"
    exit 1
}

# helper function to control the mode errors
reject_file() {
    local type="$1" f="$2" msg="$3"
    begin_file "$type" "$f"
    err "$msg"
    end_file
    return $?
}

# helper function to control fails (e.g., missing helper tool)
internal_fail_file() {
    local type="$1" f="$2" msg="$3"
    begin_file "$type" "$f"
    fail "$msg"
    end_file
    return $?
}

##########################################################################
# Options
##########################################################################


mode=""                # run or analysis
file="STDIN"                # STDIN from wrapper
samples=""             # metadata CSV file for VCF sample-name checks

#get options and arguments 
while [[ $# -gt 0 ]]; do
    case "$1" in
        -r|--run)
            mode="run"
            shift
            ;;
        -a|--analysis)
            mode="analysis"
            shift
            ;;
        -f|--file)
            [[ $# -lt 2 ]] && { echo "ERROR: --file requires an argument"; usage; }
            file="$2"
            shift 2
            ;;
        -e|--extension)
            [[ $# -lt 2 ]] && { echo "ERROR: --extension requires an argument"; usage; }
            extension="$2"
            shift 2
            ;;
        -s|--samples)
            [[ $# -lt 2 ]] && { echo "ERROR: --samples requires an argument"; usage; }
            samples="$2"
            shift 2
            ;;
        -h|--help)
            usage
            ;;
        *)
            echo "ERROR: Unknown option $1"
            usage
            ;;
    esac
done

if [[ -z "$mode" ]]; then
    echo "Error: you must specify either --run or --analysis"
    usage
fi

if [[ -z "$extension" ]]; then
    echo "Error: --extension is required"
    usage
fi

##############################################################################
# Shared status parsing
##############################################################################

print_status() {
    local level="$1"
    local message="$2"

    printf '%s\t%s\n' "$level" "$message"
}

parse_status_lines() {
    local line

    while IFS= read -r line; do
        case "$line" in
            OK$'\t'*)
                ok "${line#OK	}"
                ;;
            FAIL$'\t'*)
                fail "${line#FAIL	}"
                ;;
            ERROR$'\t'*)
                err "${line#ERROR	}"
                ;;
            "")
                ;;
            *)
                err "Unexpected QC output: $line"
                ;;
        esac
    done
}


##############################################################################
# FASTQ
##############################################################################

normalize_fastq_stream() {
    case "$extension" in
        fastq|fq)
            cat -
            ;;
        fastq.gz|fq.gz)
            gunzip -c -
            ;;
        fastq.bz2|fq.bz2)
            bunzip2 -c -
            ;;
    esac
}

fastq_basic_checks() {
    awk '
        NR == 1 {
            if ($0 !~ /^@/) {
                print "ERROR\tFormat error: record 1 does not start with @."
            }
        }

        {
            line_count++
        }

        END {
            if (line_count == 0) {
                print "ERROR\tFile is empty. Please upload a non-empty FASTQ file."
                exit 0
            }

            print "OK\tFASTQ line count: " line_count "."

            if (line_count % 4 != 0) {
                print "ERROR\tFASTQ line count is not divisible by 4. The file may be truncated or malformed."
            } else {
                print "OK\tFASTQ line count is divisible by 4."
            }
        }
    '
}

fastq_validator_check() {
    #check if fastQValidator is available
    if ! command -v fastQValidator >/dev/null 2>&1; then
        print_status "FAIL" "fastQValidator not found"
        return 0 
    fi

    local vout rc errs

    vout=$(fastQValidator --file /dev/stdin --disableSeqIDCheck 2>&1)
    rc=$?

    if (( rc == 0)); then 
        print_status "OK" "fastQValidator passed without errors."
        return 0 
    fi

    errs=$(
        printf '%s\n' "$vout" \
        | grep '^ERROR on Line ' \
        | tr '\n' '; ' \
        | sed 's/; $//'
    )

    if [[ -z "$errs" ]]; then 
        print_status "ERROR" "fastQvalidator failed, but no detailed error message was returned. Please check if FASTQ format is valid"
    else
        print_status "ERROR" "fastQvalidator failed: $errs"
    fi

    return 0

}

check_fastq_stdin() {
    local qc output rc

    begin_file "FASTQ" "$file"

    qc_output=$(
        normalize_fastq_stream \
            | tee >(fastq_validator_check) \
            | fastq_basic_checks
    )

    rc=$?

    if (( rc != 0 )); then
        err "Failed to decompress/read FASTQ file. Please check file integrity and resubmit."
        end_file
        return $?
    fi

    parse_status_lines <<< "$qc_output"

    end_file
    return $?
}


##############################################################################
# Main
##############################################################################


# -- determine file type and execute checks -- 
case "$extension" in
  fastq|fastq.gz|fastq.bz2|fq|fq.gz|fq.bz2)
    if [[ "$mode" == "run" ]]; then
        check_fastq_stdin "$file"
        exit $?
    else
        reject_file "FASTQ" "$file" "FASTQ files need to be uploaded as RUNs."
        exit $?
    fi
    ;;

#   bam|bam.gz)
#     check_bam "$file"
#     exit $?
#     ;;
#   cram|cram.gz) 
#     if [[ "$mode" == "analysis" ]]; then
#         check_cram "$file"
#         exit $?
#     else 
#         reject_file "CRAM" "$file" "CRAM files need to be uploaded as ANALYSIS"
#         exit $?
#     fi
#     ;;
    
#   vcf|vcf.gz|vcf.bz2)
#     if [[ "$mode" == "run" ]]; then
#         reject_file "VCF" "$file" "VCF/BCF files need to be uploaded as ANALYSIS"
#         exit $?
#     fi
    
#     if [[ -z "$samples" ]]; then
#         internal_fail_file "VCF" "$file" "Sample metadata CSV (-s) was not provided to the QC script"
#         exit $?
#     fi

#     check_vcf "$file"
#     exit $?
#     ;;
  *)
    echo "[WARNING] FILE $file - unsupported extension; skipping"
    exit 0
    ;;
esac