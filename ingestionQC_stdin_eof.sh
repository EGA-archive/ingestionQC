#!/usr/bin/env bash

set +e
set -o pipefail

export TERM=dumb
export NO_COLOR=1

SCRIPT_DIR="$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd -P)"

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
    echo "  -e, --extension EXT      File extension (e.g., fastq.gz, vcf.gz, bam, cram)"
    echo "  -r, --run                Run mode"
    echo "  -a, --analysis           Analysis mode"
    echo "  -f, --file STDIN         Input STDIN/display name to check"
    echo "  -s, --samples FILE       Metadata CSV file used for VCF sample-name checks"
    echo "  -h, --help               Show this help message and exit"
    echo
    echo "Environment:"
    echo "  VCF_RECORDS=N            Number of VCF records sent to VCFX_validator, default 100000"
    echo "  BAM_RECORDS=N            Number of BAM records inspected, default 100000"
    echo "  CRAM_RECORDS=N           Number of CRAM records inspected, default 100000"
    echo "  BAM_BGZF_CHECK=PATH      Path to bam_bgzf_check.py, default: same directory as this script"
    echo "  BAM_BGZF_NO_CRC=1        Faster BAM BGZF check without CRC/ISIZE validation"
    echo "  CRAM_EOF_CHECK=PATH      Path to cram_eof_check.py, default: same directory as this script"
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

print_status() {
    local level="$1"
    local message="$2"

    printf '%s\t%s\n' "$level" "$message"
}

parse_status_lines() {
    local line
    local level
    local message

    while IFS= read -r line; do
        [[ -z "$line" ]] && continue

        if [[ "$line" != *$'\t'* ]]; then
            err "Unexpected QC output: $line"
            continue
        fi

        level="${line%%$'\t'*}"
        message="${line#*$'\t'}"

        case "$level" in
            OK)
                ok "$message"
                ;;
            FAIL)
                fail "$message"
                ;;
            ERROR)
                err "$message"
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

VCF_RECORDS="${VCF_RECORDS:-100000}"
BAM_RECORDS="${BAM_RECORDS:-100000}"
CRAM_RECORDS="${CRAM_RECORDS:-100000}"

BAM_BGZF_CHECK="${BAM_BGZF_CHECK:-$SCRIPT_DIR/bam_bgzf_check.py}"
BAM_BGZF_NO_CRC="${BAM_BGZF_NO_CRC:-0}"
CRAM_EOF_CHECK="${CRAM_EOF_CHECK:-$SCRIPT_DIR/cram_eof_check.py}"

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

    validator_output=$(printf '%s\n' "$validator_output" | strip_ansi)

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
    local allowed_tmp
    local bcftools_output
    local rc

    allowed_tmp=$(mktemp) || {
        print_status "FAIL" "Could not create temporary file for metadata sample names."
        return 0
    }

    awk -F',' '
        NR > 1 {
            for (i = 1; i <= 3; i++) {
                gsub(/^[ \t\r]+|[ \t\r]+$/, "", $i)
                if ($i != "") print $i
            }
        }
    ' "$samples" | sort -u > "$allowed_tmp"

    rc=$?

    if (( rc != 0 )); then
        rm -f "$allowed_tmp"
        print_status "FAIL" "Could not read sample metadata CSV."
        return 0
    fi

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

bam_bgzf_check() {
    local bgzf_output
    local rc

    if [[ ! -x "$BAM_BGZF_CHECK" ]]; then
        print_status "FAIL" "BAM BGZF checker not found or not executable: $BAM_BGZF_CHECK."
        return 0
    fi

    if [[ "$BAM_BGZF_NO_CRC" == "1" ]]; then
        bgzf_output=$("$BAM_BGZF_CHECK" --no-crc 2>/dev/null)
    else
        bgzf_output=$("$BAM_BGZF_CHECK" 2>/dev/null)
    fi

    rc=$?

    bgzf_output=$(printf '%s\n' "$bgzf_output" | strip_ansi)

    if [[ -n "$bgzf_output" ]]; then
        printf '%s\n' "$bgzf_output"
        return 0
    fi

    if (( rc != 0 )); then
        print_status "ERROR" "BAM BGZF/EOF check failed without a detailed error message."
    else
        print_status "OK" "BAM BGZF/EOF check passed."
    fi

    return 0
}

bam_samtools_check() {
    local check_output
    local rc

    check_output=$(
        samtools view -h - 2>&1 \
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
                    samtools_error = ""
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

                /^[[]/ || /^samtools/ || /^E::/ || /^W::/ {
                    samtools_error = samtools_error $0 "; "
                    next
                }

                {
                    records_seen++

                    flag = $2 + 0

                    if (!has_flag(flag, 4) && !has_flag(flag, 256) && !has_flag(flag, 2048)) {
                        primary_mapped_seen = 1
                    }

                    if (records_seen >= max_records) {
                        exit
                    }
                }

                END {
                    if (samtools_error != "") {
                        sub(/; $/, "", samtools_error)
                        print "ERROR\tsamtools reported an error while reading BAM stream: " samtools_error
                        exit 0
                    }

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
                        print "OK\tBAM appears to contain primary mapped alignments in the first " records_seen " inspected records."

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
                        print "OK\tNo primary mapped alignments detected in the first " records_seen " inspected records."

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

    if (( rc != 0 )) && [[ -z "$check_output" ]]; then
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

    if [[ ! -x "$BAM_BGZF_CHECK" ]]; then
        fail "BAM BGZF checker not found or not executable: $BAM_BGZF_CHECK."
        end_file
        return $?
    fi

    qc_output=$(
        {
            cat - \
                | tee -p \
                    >(bam_bgzf_check >&3) \
                    >(bam_refgen_check >&3) \
                    >(bam_samtools_check >&3) \
                > /dev/null
        } 3>&1
    )
    rc=$?

    qc_output=$(printf '%s\n' "$qc_output" | strip_ansi)

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
# CRAM
##############################################################################

cram_eof_check() {
    local eof_output
    local rc

    if [[ ! -x "$CRAM_EOF_CHECK" ]]; then
        print_status "FAIL" "CRAM EOF checker not found or not executable: $CRAM_EOF_CHECK."
        return 0
    fi

    eof_output=$("$CRAM_EOF_CHECK" 2>/dev/null)
    rc=$?

    eof_output=$(printf '%s\n' "$eof_output" | strip_ansi)

    if [[ -n "$eof_output" ]]; then
        printf '%s\n' "$eof_output"
        return 0
    fi

    if (( rc != 0 )); then
        print_status "ERROR" "CRAM EOF check failed without a detailed error message."
    else
        print_status "OK" "CRAM EOF check passed."
    fi

    return 0
}

cram_samtools_check() {
    local check_output
    local rc

    check_output=$(
        samtools view -h - 2>&1 \
            | awk -v mode="$mode" -v max_records="$CRAM_RECORDS" '
                function has_flag(flag, bit) {
                    return int(flag / bit) % 2
                }

                BEGIN {
                    header_seen = 0
                    sq_seen = 0
                    sq_without_m5 = 0
                    so = ""
                    records_seen = 0
                    primary_mapped_seen = 0
                    samtools_error = ""
                }

                /^@/ {
                    header_seen = 1

                    if ($1 == "@SQ") {
                        sq_seen = 1

                        if ($0 !~ /(^|[ \t])M5:/) {
                            sq_without_m5++
                        }
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

                /^[[]/ || /^samtools/ || /^E::/ || /^W::/ {
                    samtools_error = samtools_error $0 "; "
                    next
                }

                {
                    records_seen++

                    flag = $2 + 0

                    if (!has_flag(flag, 4) && !has_flag(flag, 256) && !has_flag(flag, 2048)) {
                        primary_mapped_seen = 1
                    }

                    if (records_seen >= max_records) {
                        exit
                    }
                }

                END {
                    if (samtools_error != "") {
                        sub(/; $/, "", samtools_error)
                        print "ERROR\tsamtools reported an error while reading CRAM stream: " samtools_error
                        exit 0
                    }

                    if (!header_seen) {
                        print "ERROR\tCRAM header is missing or unreadable."
                        exit 0
                    }

                    if (!sq_seen) {
                        print "ERROR\tCRAM header is missing @SQ reference sequence entries."
                    } else {
                        print "OK\tCRAM header contains @SQ reference sequence entries."
                    }

                    if (sq_without_m5 > 0) {
                        print "ERROR\tCRAM header has @SQ reference sequence entries without M5 reference MD5 tags."
                    } else if (sq_seen) {
                        print "OK\tAll CRAM @SQ reference sequence entries contain M5 tags."
                    }

                    if (mode != "analysis") {
                        print "ERROR\tCRAM files need to be uploaded as ANALYSIS."
                    } else {
                        print "OK\tCRAM uploaded as ANALYSIS."
                    }

                    if (primary_mapped_seen) {
                        print "OK\tCRAM appears to contain primary mapped alignments in the first " records_seen " inspected records."

                        if (so == "coordinate") {
                            print "OK\tAligned CRAM header reports coordinate sorting."
                        } else {
                            print "ERROR\tAligned CRAM is not marked as coordinate sorted in the header (SO:" (so == "" ? "missing" : so) ")."
                        }
                    } else {
                        print "ERROR\tNo primary mapped alignments detected in the first " records_seen " inspected records. CRAM files are expected to be uploaded as ANALYSIS."
                    }
                }
            '
    )
    rc=$?

    if (( rc != 0 )) && [[ -z "$check_output" ]]; then
        print_status "ERROR" "samtools failed to read CRAM stream. Please check file integrity and CRAM format."
        return 0
    fi

    printf '%s\n' "$check_output"
    return 0
}

cram_refgen_check() {
    local rfg_output
    local rc
    local species
    local reference

    rfg_output=$(refgenDetector -f - -t BAM/CRAM 2>&1)
    rc=$?

    rfg_output=$(printf '%s\n' "$rfg_output" | strip_ansi)

    if (( rc != 0 )); then
        print_status "ERROR" "refgenDetector failed to inspect CRAM stream."
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

check_cram_stdin() {
    local qc_output
    local rc

    begin_file "CRAM" "$file"

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

    if [[ ! -x "$CRAM_EOF_CHECK" ]]; then
        fail "CRAM EOF checker not found or not executable: $CRAM_EOF_CHECK."
        end_file
        return $?
    fi

    qc_output=$(
        {
            cat - \
                | tee -p \
                    >(cram_eof_check >&3) \
                    >(cram_refgen_check >&3) \
                    >(cram_samtools_check >&3) \
                > /dev/null
        } 3>&1
    )
    rc=$?

    qc_output=$(printf '%s\n' "$qc_output" | strip_ansi)

    if (( rc != 0 )) && [[ -z "$qc_output" ]]; then
        err "CRAM stream processing failed before QC results could be collected. Please check file integrity and CRAM format."
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

    bam)
        check_bam_stdin
        exit $?
        ;;

    cram)
        check_cram_stdin
        exit $?
        ;;

    *)
        echo "[WARNING] FILE $file - unsupported extension; skipping"
        exit 0
        ;;
esac