#!/usr/bin/env python3

import argparse
import struct
import sys
import time
import zlib


BAM_BGZF_EOF = bytes.fromhex(
    "1f8b08040000000000ff0600424302001b0003000000000000000000"
)


def read_exact(stream, n):
    chunks = []
    total = 0

    while total < n:
        chunk = stream.read(n - total)
        if not chunk:
            break

        chunks.append(chunk)
        total += len(chunk)

    return b"".join(chunks)


def parse_bgzf_bsize(extra):
    pos = 0

    while pos + 4 <= len(extra):
        subfield_id = extra[pos:pos + 2]
        subfield_len = struct.unpack("<H", extra[pos + 2:pos + 4])[0]
        pos += 4

        if pos + subfield_len > len(extra):
            raise ValueError("malformed BGZF extra field")

        payload = extra[pos:pos + subfield_len]
        pos += subfield_len

        if subfield_id == b"BC":
            if subfield_len != 2:
                raise ValueError("malformed BGZF BC subfield")

            return struct.unpack("<H", payload)[0] + 1

    raise ValueError("missing BGZF BC subfield")


def verify_deflate_payload(block, xlen):
    payload_start = 12 + xlen
    payload_end = len(block) - 8

    if payload_end < payload_start:
        raise ValueError("BGZF block is too short to contain a gzip trailer")

    compressed_payload = block[payload_start:payload_end]
    trailer = block[-8:]

    expected_crc, expected_isize = struct.unpack("<II", trailer)

    try:
        uncompressed = zlib.decompress(compressed_payload, wbits=-15)
    except zlib.error as exc:
        raise ValueError(f"BGZF deflate payload could not be decompressed: {exc}") from exc

    observed_crc = zlib.crc32(uncompressed) & 0xFFFFFFFF
    observed_isize = len(uncompressed) & 0xFFFFFFFF

    if observed_crc != expected_crc:
        raise ValueError("BGZF block CRC check failed")

    if observed_isize != expected_isize:
        raise ValueError("BGZF block uncompressed-size check failed")


def check_bgzf_stream(stream, verify_crc=True):
    offset = 0
    block_count = 0
    last_block = b""
    seen_eof_block = False

    while True:
        header = read_exact(stream, 12)

        if header == b"":
            break

        if len(header) < 12:
            return False, (
                f"BGZF block header is truncated at byte offset {offset}; "
                f"expected 12 bytes, found {len(header)}."
            )

        if seen_eof_block:
            return False, (
                f"Additional data found after BAM EOF marker at byte offset {offset}."
            )

        if header[0:2] != b"\x1f\x8b":
            return False, f"Invalid gzip/BGZF magic bytes at byte offset {offset}."

        if header[2] != 8:
            return False, f"Invalid gzip compression method at byte offset {offset}."

        if not (header[3] & 0x04):
            return False, (
                f"BGZF block at byte offset {offset} is missing the FEXTRA flag."
            )

        if header[3] != 0x04:
            return False, (
                f"Unsupported gzip flags in BGZF block at byte offset {offset}: "
                f"0x{header[3]:02x}."
            )

        xlen = struct.unpack("<H", header[10:12])[0]
        extra = read_exact(stream, xlen)

        if len(extra) < xlen:
            return False, (
                f"BGZF extra field is truncated at byte offset {offset}; "
                f"expected {xlen} bytes, found {len(extra)}."
            )

        try:
            block_size = parse_bgzf_bsize(extra)
        except ValueError as exc:
            return False, f"Malformed BGZF block at byte offset {offset}: {exc}."

        already_read = 12 + xlen
        remaining = block_size - already_read

        if remaining < 8:
            return False, (
                f"Malformed BGZF block at byte offset {offset}: "
                f"block size {block_size} is too small."
            )

        rest = read_exact(stream, remaining)

        if len(rest) < remaining:
            return False, (
                f"BGZF block is truncated at byte offset {offset}; "
                f"expected block size {block_size} bytes, "
                f"found {already_read + len(rest)} bytes."
            )

        block = header + extra + rest

        if verify_crc:
            try:
                verify_deflate_payload(block, xlen)
            except ValueError as exc:
                return False, f"Invalid BGZF block at byte offset {offset}: {exc}."

        block_count += 1
        last_block = block

        if block == BAM_BGZF_EOF:
            seen_eof_block = True

        offset += block_size

    if block_count == 0:
        return False, "BAM stream is empty. Please upload a non-empty BAM file."

    if last_block != BAM_BGZF_EOF:
        return False, (
            "BAM/BGZF EOF marker is missing or invalid. "
            "The stream may be incomplete or truncated."
        )

    return True, (
        f"BAM/BGZF stream structure is valid; EOF marker is present; "
        f"BGZF blocks checked={block_count}."
    )


def emit_status(level, message):
    sys.stdout.write(f"{level}\t{message}\n")


def main():
    parser = argparse.ArgumentParser(
        description="Check BAM/BGZF block structure and BAM EOF marker from STDIN."
    )
    parser.add_argument(
        "--no-crc",
        action="store_true",
        help=(
            "Only check BGZF block structure and final BAM EOF marker. "
            "Do not decompress blocks to verify CRC/ISIZE."
        ),
    )
    parser.add_argument(
        "--show-runtime",
        action="store_true",
        help="Print runtime to stderr.",
    )

    args = parser.parse_args()

    start = time.monotonic()

    ok, message = check_bgzf_stream(
        sys.stdin.buffer,
        verify_crc=not args.no_crc,
    )

    elapsed = time.monotonic() - start

    if args.show_runtime:
        sys.stderr.write(f"Runtime: {elapsed:.3f} seconds\n")

    if ok:
        emit_status("OK", message)
        return 0

    emit_status("ERROR", message)
    return 1


if __name__ == "__main__":
    sys.exit(main())