#!/usr/bin/env bash

set +e
set -o pipefail

# Force plain-text tool output where possible.
export TERM=dumb
export NO_COLOR=1

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

ok()   { _OKS+=("$1"); }
fail() { _FAILS+=("$1"); }
err()  { _ERRS+=("$1"); }

end_file() {
    local msg

    if ((${#_ERRS[@]})); then
        msg=$(printf '%s; ' "${_ERRS[@]}")
        msg=${msg%; }
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
    echo "Usage: $0 [OPTIONS] -f <stdin> -e <extension> [-s <samples>]"
    echo "Options:"
    echo "  -e, --extension EXT      File extension (e.g., fastq.gz, vcf.gz)"
    echo "  -r, --run                Run mode"
    echo "  -a, --analysis           Analysis mode"
    echo "  -f, --file STDIN         Input STDIN/display name to check"
    echo "  -s, --samples FILE       Metadata CSV file used for VCF sample-name checks"
    echo "  --debug                  Print debug messages to stderr"
    echo "  -h, --help               Show this help message and exit"
    echo
    echo "Environment:"
    echo "  VCF_RECORDS=N            Number of VCF records sent to VCFX_validator, default 10000"
    exit 1
}

reject_file() {
    local type="$1" f="$2" msg="$3"
    begin_file "$type" "$f"
    err "$msg"
    end_file
    return $?
}

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

debug_log() {
    if [[ "$debug" -eq 1 ]]; then
        printf '[DEBUG] %s\n' "$*" >&2
    fi
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

##############################################################################
# Options
##############################################################################

extension=""
mode=""
file="STDIN"
samples=""
debug=0
VCF_RECORDS="${VCF_RECORDS:-10000}"

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
        --debug)
            debug=1
            shift
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

    debug_log "check_fastq_stdin: started extension=$extension mode=$mode file=$file"

    if ! command -v fastQValidator >/dev/null 2>&1; then
        fail "fastQValidator not found."
        debug_log "check_fastq_stdin: fastQValidator not found"
        end_file
        return $?
    fi

    debug_log "check_fastq_stdin: running fastQValidator"

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

    validator_output=$(printf '%s\n' "$validator_output" | strip_ansi)

    debug_log "check_fastq_stdin: fastQValidator pipeline finished rc=$rc"

    if [[ "$debug" -eq 1 ]]; then
        printf '%s\n' "$validator_output" | sed 's/^/[DEBUG] FASTQ_VALIDATOR_OUTPUT: /' >&2
    fi

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

    ok "fastQValidator passed."
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
    debug_log "vcf_sample_check: started"

    local allowed_tmp
    local bcftools_output
    local rc
    local allowed_count

    allowed_tmp=$(mktemp) || {
        print_status "FAIL" "Could not create temporary file for metadata sample names."
        debug_log "vcf_sample_check: failed to create temporary file"
        return 0
    }

    debug_log "vcf_sample_check: created temp file $allowed_tmp"

    awk -F',' '
        NR > 1 {
            for (i = 1; i <= 3; i++) {
                gsub(/^[ \t\r]+|[ \t\r]+$/, "", $i)
                if ($i != "") print $i
            }
        }
    ' "$samples" | sort -u > "$allowed_tmp"

    rc=$?
    allowed_count=$(wc -l < "$allowed_tmp" 2>/dev/null)

    debug_log "vcf_sample_check: loaded ${allowed_count:-0} allowed sample identifiers from $samples"

    if (( rc != 0 )); then
        rm -f "$allowed_tmp"
        print_status "FAIL" "Could not read sample metadata CSV."
        debug_log "vcf_sample_check: metadata CSV parsing failed rc=$rc"
        return 0
    fi

    if [[ ! -s "$allowed_tmp" ]]; then
        rm -f "$allowed_tmp"
        print_status "FAIL" "No sample identifiers could be read from metadata CSV."
        debug_log "vcf_sample_check: allowed sample list is empty"
        return 0
    fi

    debug_log "vcf_sample_check: running bcftools query -l"

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
                    if ($0 != "") {
                        vcf_samples++
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

    debug_log "vcf_sample_check: bcftools query branch finished rc=$rc"

    rm -f "$allowed_tmp"
    debug_log "vcf_sample_check: removed temp file"

    if (( rc != 0 )); then
        print_status "ERROR" "Could not extract sample names from VCF stream with bcftools."
        return 0
    fi

    printf '%s\n' "$bcftools_output"

    debug_log "vcf_sample_check: finished"
    return 0
}

vcf_validator_check() {
    debug_log "vcf_validator_check: started"

    local vout
    local rc
    local errs

    if ! command -v VCFX_validator >/dev/null 2>&1; then
        print_status "FAIL" "VCFX_validator not found."
        debug_log "vcf_validator_check: VCFX_validator not found"
        return 0
    fi

    debug_log "vcf_validator_check: validating header plus first $VCF_RECORDS records"

    vout=$(vcf_prefix_for_validator | VCFX_validator 2>&1)
    rc=$?

    debug_log "vcf_validator_check: VCFX_validator finished rc=$rc"

    vout=$(printf '%s\n' "$vout" | strip_ansi)

    if (( rc == 0 )) || printf '%s\n' "$vout" | grep -q '^Status:[[:space:]]*PASSED'; then
        print_status "OK" "VCF passed VCFX_validator on header plus first ${VCF_RECORDS} records."
        debug_log "vcf_validator_check: passed"
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

    debug_log "vcf_validator_check: finished with validation error"
    return 0
}

check_vcf_stdin() {
    local qc_output
    local rc

    begin_file "VCF" "$file"

    debug_log "check_vcf_stdin: started extension=$extension mode=$mode file=$file samples=$samples VCF_RECORDS=$VCF_RECORDS"

    if [[ "$mode" != "analysis" ]]; then
        err "VCF files need to be uploaded as ANALYSIS."
        debug_log "check_vcf_stdin: rejected mode=$mode"
        end_file
        return $?
    fi

    if [[ -z "$samples" ]]; then
        fail "Sample metadata CSV (-s) was not provided to the QC script."
        debug_log "check_vcf_stdin: samples argument missing"
        end_file
        return $?
    fi

    if [[ ! -f "$samples" ]]; then
        fail "Sample metadata CSV '$samples' was not found."
        debug_log "check_vcf_stdin: samples file not found: $samples"
        end_file
        return $?
    fi

    if ! command -v bcftools >/dev/null 2>&1; then
        fail "bcftools not found; VCF sample names cannot be checked."
        debug_log "check_vcf_stdin: bcftools not found"
        end_file
        return $?
    fi

    if ! command -v VCFX_validator >/dev/null 2>&1; then
        fail "VCFX_validator not found."
        debug_log "check_vcf_stdin: VCFX_validator not found"
        end_file
        return $?
    fi

    debug_log "check_vcf_stdin: starting tee pipeline"

    qc_output=$(
        {
            normalize_vcf_stream \
                | tee -p >(vcf_sample_check >&3) \
                | vcf_validator_check
        } 3>&1
    )
    rc=$?

    debug_log "check_vcf_stdin: tee pipeline finished rc=$rc"

    qc_output=$(printf '%s\n' "$qc_output" | strip_ansi)

    if [[ "$debug" -eq 1 ]]; then
        debug_log "check_vcf_stdin: collected QC output follows"
        printf '%s\n' "$qc_output" | sed 's/^/[DEBUG] QC_OUTPUT: /' >&2
    fi

    if (( rc != 0 )); then
        err "VCF stream processing failed before QC results could be collected. Please check file integrity and VCF format."
        debug_log "check_vcf_stdin: failing because pipeline rc=$rc"
        end_file
        return $?
    fi

    parse_status_lines <<< "$qc_output"

    debug_log "check_vcf_stdin: parsed status lines"

    end_file
    return $?
}

##############################################################################
# Main
##############################################################################

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

    vcf|vcf.gz|vcf.bz2)
        if [[ "$mode" == "run" ]]; then
            reject_file "VCF" "$file" "VCF files need to be uploaded as ANALYSIS."
            exit $?
        fi

        check_vcf_stdin
        exit $?
        ;;

    *)
        echo "[WARNING] FILE $file - unsupported extension; skipping"
        exit 0
        ;;
esac