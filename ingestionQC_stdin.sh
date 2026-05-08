#!/usr/bin/env bash

set +e  # aggregates errors; we do not want to abort afer the first failure but collect all possible issues
set -o pipefail # ensures that if any command in a pipeline fails, the entire pipeline is considered to have failed (i.e., it will return a non-zero exit status). This is important for error handling, especially when using tools like grep or awk in pipelines, as it allows us to detect failures that might otherwise be masked by successful commands later in the pipeline.

#forces a non-interactive terminal mode and disables colored output, which is important for consistent parsing of error messages from tools like fastQValidator, especially when stripping ANSI escape codes. This ensures that the output is plain text without any formatting characters that could interfere with error detection and reporting.
export TERM=dumb
export NO_COLOR=1

declare -a _OKS _FAILS _ERRS
_cur_type=""
_cur_file=""

VCF_RECORDS=100000

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

##############################################################################
# Shared helpers
##############################################################################

strip_ansi() {
    sed -E $'s/\x1B\\[[0-9;?]*[A-Za-z]//g'
}

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


##########################################################################
# Options
##########################################################################

extension=""         # file extension (fastq.gz, bam, cram, vcf.gz, etc.)
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

strip_ansi() {
    sed -E $'s/\x1B\\[[0-9;?]*[A-Za-z]//g'
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

check_fastq_stdin() {
    local validator_output
    local rc
    local errs

    begin_file "FASTQ" "$file"

    if ! command -v fastQValidator >/dev/null 2>&1; then
        fail "fastQValidator not found."
        end_file
        return $?
    fi

    validator_output=$(
        {
            normalize_fastq_stream \
                | awk '
                    {
                        seen = 1
                        print
                    }

                    END {
                        if (!seen) {
                            print "__INGESTIONQC_EMPTY_FASTQ__" > "/dev/stderr"
                            exit 1
                        }
                    }
                ' \
                | fastQValidator --file - --disableSeqIDCheck
        } 2>&1
    )
    rc=$?

    if (( rc != 0 )); then
        if printf '%s\n' "$validator_output" | grep -q '__INGESTIONQC_EMPTY_FASTQ__'; then
            err "File is empty. Please upload a non-empty FASTQ file."
        else
            errs=$(
                printf '%s\n' "$validator_output" \
                    | tr -d '\r' \
                    | grep -E '^ERROR' \
                    | paste -sd ';' -
            )

            if [[ -z "$errs" ]]; then
                err "Failed to decompress/read FASTQ file or fastQValidator failed without a detailed error message. Please check file integrity and FASTQ format."
            else
                err "fastQValidator failed: $errs"
            fi
        fi

        end_file
        return $?
    fi

    ok "fastQValidator passed and FASTQ line count is valid."
    end_file
    return $?
}


##############################################################################
# VCF
##############################################################################

normalize_vcf_stream() {
    case "$extension" in
        vcf)
            cat -
            ;;
        vcf.gz)
            gunzip -c -
            ;;
        vcf.bz2)
            bunzip2 -c -
            ;;
    esac
}

vcf_prefix_for_validator() {
    awk -v max_records="$VCF_RECORDS" '
        /^#/ {
            print
            next
        }

        records < max_records {
            print
            records++
            next
        }

        records >= max_records {
            exit
        }
    '
}

vcf_sample_check() {
    local allowed_tmp
    local bcftools_output
    local rc

    allowed_tmp=$(mktemp) || {
        print_status "FAIL" "Could not create temporary file for metadata sample names."
        return 0
    }

    # Metadata CSV columns:
    #   alias,stable_id,sample_id
    #
    # A VCF sample is accepted if it matches any non-empty value
    # in alias, stable_id, or sample_id.
    awk -F',' '
        NR > 1 {
            for (i = 1; i <= 3; i++) {
                gsub(/^[ \t\r]+|[ \t\r]+$/, "", $i)
                if ($i != "") print $i
            }
        }
    ' "$samples" | sort -u > "$allowed_tmp"

    if [[ ! -s "$allowed_tmp" ]]; then
        rm -f "$allowed_tmp"
        print_status "FAIL" "No sample identifiers could be read from metadata CSV."
        return 0
    fi

    bcftools_output=$(
        bcftools query -l 2>&1 \
            | awk -v allowed_file="$allowed_tmp" '
                BEGIN {
                    while ((getline s < allowed_file) > 0) {
                        allowed[s] = 1
                    }
                    close(allowed_file)
                }

                {
                    if ($0 != "" && !($0 in allowed)) {
                        missing[$0] = 1
                    }
                }

                END {
                    count = 0
                    list = ""

                    for (s in missing) {
                        count++
                        if (count <= 20) {
                            list = list (list == "" ? "" : ",") s
                        }
                    }

                    if (count > 0) {
                        if (count > 20) {
                            print "ERROR\tSamples present in VCF but missing from registered metadata: " list ", ... (" count " total missing samples)."
                        } else {
                            print "ERROR\tSamples present in VCF but missing from registered metadata: " list "."
                        }
                    } else {
                        print "OK\tAll VCF samples are present in registered metadata."
                    }
                }
            '
    )
    rc=$?

    rm -f "$allowed_tmp"

    if (( rc != 0 )); then
        print_status "ERROR" "Could not extract sample names from VCF stream with bcftools."
        return 0
    fi

    printf '%s\n' "$bcftools_output"
    return 0
}


vcf_validator_check() {
    local vout
    local rc
    local errs

    if ! command -v VCFX_validator >/dev/null 2>&1; then
        print_status "FAIL" "VCFX_validator not found."
        return 0
    fi

    vout=$(vcf_prefix_for_validator | VCFX_validator 2>&1)
    rc=$?

    vout=$(printf '%s\n' "$vout" | strip_ansi)

    if (( rc == 0 )) || printf '%s\n' "$vout" | grep -q '^Status:[[:space:]]*PASSED'; then
        print_status "OK" "VCF passed VCFX_validator on header plus first ${VCF_RECORDS} records."
        return 0
    fi

    errs=$(
        printf '%s\n' "$vout" \
            | tr -d '\r' \
            | grep -E '^(ERROR|Error|error):?' \
            | paste -sd ';' -
    )

    if [[ -z "$errs" ]]; then
        errs=$(
            printf '%s\n' "$vout" \
                | grep -v '^Status:[[:space:]]*PASSED' \
                | grep -v '^Lines read' \
                | grep -v '^$' \
                | head -n 10 \
                | paste -sd ';' -
        )
    fi

    if [[ -z "$errs" ]]; then
        print_status "ERROR" "VCFX_validator failed on header plus first ${VCF_RECORDS} records, but no detailed error message was returned. Please check VCF format."
    else
        print_status "ERROR" "VCF validation failed on header plus first ${VCF_RECORDS} records: $errs"
    fi

    return 0
}


check_vcf_stdin() {
    local qc_output
    local rc

    begin_file "VCF" "$file"

    if [[ "$mode" != "analysis" ]]; then
        err "VCF files need to be uploaded as ANALYSIS."
        end_file
        return $?
    fi

    if [[ -z "$samples" ]]; then
        fail "Sample metadata CSV (-s) was not provided to the QC script."
        end_file
        return $?
    fi

    if [[ ! -f "$samples" ]]; then
        fail "Sample metadata CSV '$samples' was not found."
        end_file
        return $?
    fi

    if ! command -v bcftools >/dev/null 2>&1; then
        fail "bcftools not found; VCF sample names cannot be checked."
        end_file
        return $?
    fi

    if ! command -v VCFX_validator >/dev/null 2>&1; then
        fail "VCFX_validator not found."
        end_file
        return $?
    fi

    qc_output=$(
        {
            normalize_vcf_stream \
                | tee -p >(vcf_sample_check >&3) \
                | vcf_validator_check
        } 3>&1
    )
    rc=$?

    qc_output=$(printf '%s\n' "$qc_output" | strip_ansi)

    # rc=141 is expected when tee/gunzip receives SIGPIPE because one or more
    # QC branches finished early. This is not fatal if we collected QC output.
    if (( rc != 0 )) && [[ -z "$qc_output" ]]; then
        err "VCF stream processing failed before QC results could be collected. Please check file integrity and VCF format."
        end_file
        return $?
    fi

    parse_status_lines <<< "$qc_output"

    end_file
    return $?
}

##############################################################################
# BAM
##############################################################################

BAM_RECORDS="${BAM_RECORDS:-100000}"

bam_samtools_check() {
    local check_output
    local rc

    check_output=$(
        samtools view -h - 2>/dev/null \
            | awk -v mode="$mode" -v max_records="$BAM_RECORDS" '
                function has_flag(flag, bit) {
                    return int(flag / bit) % 2
                }

                BEGIN {
                    header_seen = 0
                    sq_seen = 0
                    so = ""
                    records_seen = 0
                    primary_mapped_seen = 0
                }

                /^@/ {
                    header_seen = 1

                    if ($1 == "@SQ") {
                        sq_seen = 1
                    }

                    if ($1 == "@HD") {
                        for (i = 1; i <= NF; i++) {
                            if ($i ~ /^SO:/) {
                                so = substr($i, 4)
                            }
                        }
                    }

                    next
                }

                {
                    records_seen++

                    flag = $2 + 0

                    # Primary mapped alignment:
                    #   not unmapped       0x4
                    #   not secondary      0x100
                    #   not supplementary  0x800
                    if (!has_flag(flag, 4) && !has_flag(flag, 256) && !has_flag(flag, 2048)) {
                        primary_mapped_seen = 1
                    }

                    if (records_seen >= max_records) {
                        exit
                    }
                }

                END {
                    if (!header_seen) {
                        print "ERROR\tBAM header is missing or unreadable."
                        exit 0
                    }

                    if (!sq_seen) {
                        print "ERROR\tBAM header is missing @SQ reference sequence entries."
                    } else {
                        print "OK\tBAM header contains @SQ reference sequence entries."
                    }

                    if (primary_mapped_seen) {
                        print "OK\tBAM appears to contain primary mapped alignments in the inspected records."

                        if (mode == "run") {
                            print "ERROR\tThis BAM appears to be ALIGNED; please upload this file as an ANALYSIS, not a RUN."
                        } else {
                            print "OK\tAligned BAM uploaded as ANALYSIS."
                        }

                        if (so == "coordinate") {
                            print "OK\tAligned BAM header reports coordinate sorting."
                        } else {
                            print "ERROR\tAligned BAM is not marked as coordinate sorted in the header (SO:" (so == "" ? "missing" : so) ")."
                        }
                    } else {
                        print "OK\tNo primary mapped alignments detected in the first " max_records " inspected records."

                        if (mode == "run") {
                            print "OK\tUnaligned BAM uploaded as RUN."
                        } else {
                            print "ERROR\tThis BAM appears to be UNALIGNED; please upload this file as a RUN, not an ANALYSIS."
                        }
                    }
                }
            '
    )
    rc=$?

    if (( rc != 0 )); then
        print_status "ERROR" "samtools failed to read BAM stream. Please check file integrity and BAM format."
        return 0
    fi

    printf '%s\n' "$check_output"
    return 0
}

bam_refgen_check() {
    local rfg_output
    local rc
    local species
    local reference

    rfg_output=$(refgenDetector -f - -t BAM/CRAM 2>&1)
    rc=$?

    rfg_output=$(printf '%s\n' "$rfg_output" | strip_ansi)

    if (( rc != 0 )); then
        print_status "ERROR" "refgenDetector failed to inspect BAM stream."
        return 0
    fi

    species=$(
        printf '%s\n' "$rfg_output" \
            | awk -F'Species detected:[[:space:]]*' '/Species detected:/ {print $2; exit}' \
            | xargs
    )

    reference=$(
        printf '%s\n' "$rfg_output" \
            | awk -F'Reference genome version[[:space:]]*:[[:space:]]*' '/Reference genome version/ {print $2; exit}' \
            | xargs
    )

    if [[ -z "$species" ]]; then
        print_status "ERROR" "refgenDetector produced no species result."
        return 0
    fi

    if [[ "$species" != "Homo sapiens" ]]; then
        print_status "ERROR" "refgenDetector: species is not human ($species)."
        return 0
    fi

    if [[ -n "$reference" ]]; then
        print_status "OK" "refgenDetector detected Homo sapiens reference genome ($reference)."
    else
        print_status "OK" "refgenDetector detected Homo sapiens."
    fi

    return 0
}

check_bam_stdin() {
    local qc_output
    local rc

    begin_file "BAM" "$file"

    if ! command -v samtools >/dev/null 2>&1; then
        fail "samtools not found."
        end_file
        return $?
    fi

    if ! command -v refgenDetector >/dev/null 2>&1; then
        fail "refgenDetector not found."
        end_file
        return $?
    fi

    qc_output=$(
        {
            cat - \
                | tee -p >(bam_refgen_check >&3) \
                | bam_samtools_check
        } 3>&1
    )
    rc=$?

    qc_output=$(printf '%s\n' "$qc_output" | strip_ansi)

    # rc=141 is expected when tee receives SIGPIPE because one QC branch
    # finished earlier than the other. This is not fatal if QC output exists.
    if (( rc != 0 )) && [[ -z "$qc_output" ]]; then
        err "BAM stream processing failed before QC results could be collected. Please check file integrity and BAM format."
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
    vcf|vcf.gz|vcf.bz2|bcf)
    if [[ "$mode" == "run" ]]; then
        reject_file "VCF" "$file" "VCF/BCF files need to be uploaded as ANALYSIS."
        exit $?
    fi

    check_vcf_stdin
    exit $?
    ;;
    bam)
        check_bam_stdin
        exit $?
        ;;
  *)
    echo "[WARNING] FILE $file - unsupported extension; skipping"
    exit 0
    ;;
esac





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