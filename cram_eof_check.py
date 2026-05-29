#!/usr/bin/env python3

import argparse
import sys
import time


CRAM_EOF_CONTAINER = bytes.fromhex(
    "0f000000ffffffff0fe0454f4600000000010005bdd94f0001000606010001000100ee63014b"
)

CHUNK_SIZE = 1024 * 1024


def emit_status(level: str, message: str) -> None:
    sys.stdout.write(f"{level}\t{message}\n")


def main() -> int:
    parser = argparse.ArgumentParser(
        description="Check CRAM EOF container from STDIN."
    )
    parser.add_argument(
        "--show-runtime",
        action="store_true",
        help="Print runtime to stderr.",
    )

    args = parser.parse_args()

    start = time.monotonic()

    tail = b""
    seen_data = False

    while True:
        chunk = sys.stdin.buffer.read(CHUNK_SIZE)

        if not chunk:
            break

        seen_data = True
        tail = (tail + chunk)[-len(CRAM_EOF_CONTAINER):]

    elapsed = time.monotonic() - start

    if args.show_runtime:
        sys.stderr.write(f"Runtime: {elapsed:.3f} seconds\n")

    if not seen_data:
        emit_status("ERROR", "CRAM stream is empty. Please upload a non-empty CRAM file.")
        return 1

    if tail != CRAM_EOF_CONTAINER:
        emit_status(
            "ERROR",
            "CRAM EOF container is missing or invalid. The stream may be incomplete or truncated.",
        )
        return 1

    emit_status("OK", "CRAM EOF container is present.")
    return 0


if __name__ == "__main__":
    sys.exit(main())