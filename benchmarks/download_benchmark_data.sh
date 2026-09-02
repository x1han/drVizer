#!/usr/bin/env bash
# download_benchmark_data.sh -- fetch inputs used by benchmarks/*.py
#
# Most drVizer benchmarks generate their own synthetic fixtures via
# the ``_write_*`` helpers in ``bench_bed_concurrency.py`` /
# ``bench_cache_vs_no_cache.py``. This script exists for the
# optional public-dataset workflows (samtools-based BAM slicing and
# real GTF / BED inputs from third parties) that a maintainer might
# run locally.
#
# Optional dependencies
# ---------------------
# This script is INTENTIONALLY graceful when optional tools are
# missing. A maintainer running it on a workstation without
# ``samtools`` / ``curl`` / ``wget`` should see clear "skipping"
# messages, not a hard failure.
#
#   * ``samtools`` -- required only for the BAM slicing step. If
#     ``samtools`` is not on ``$PATH``, the BAM step prints a skip
#     notice and exits 0. Most benchmarks do not require BAMs, so
#     this is rarely needed.
#   * ``curl`` OR ``wget`` -- required for the GTF / BED fetch
#     steps. We probe both and use whichever is present; if neither
#     is available the public-dataset steps skip with a notice.
#
# Usage
# -----
#   ./benchmarks/download_benchmark_data.sh [DEST_DIR]
#
# Defaults to ``./benchmark_data/`` under the repository root.

set -euo pipefail

REPO_ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
DEST_DIR="${1:-${REPO_ROOT}/benchmark_data}"

mkdir -p "${DEST_DIR}"

log()  { printf '[bench-data] %s\n' "$*"; }
skip() { log "skipping: $*"; }

# Pick curl or wget; fall back to a clear skip notice.
if command -v curl >/dev/null 2>&1; then
    fetch() {
        curl --fail --silent --show-error --location --output "$2" "$1"
    }
elif command -v wget >/dev/null 2>&1; then
    fetch() {
        wget --quiet --output-document="$2" "$1"
    }
else
    fetch() {
        skip "neither curl nor wget on PATH; cannot fetch $1"
        return 1
    }
fi

# 1. BAM slicing via samtools (optional).
if [[ "${SKIP_BAM:-0}" == "1" ]]; then
    skip "SKIP_BAM=1 set; skipping BAM fetch"
elif command -v samtools >/dev/null 2>&1; then
    log "samtools detected; BAM fetch path available"
    log "(no public BAM URL is hardcoded -- wire your own in here)"
else
    skip "samtools not on PATH; skipping BAM step (most benchmarks do not need it)"
fi

# 2. Synthetic GTF / BED fixtures (no network required).
SYNTH_DIR="${DEST_DIR}/synthetic"
mkdir -p "${SYNTH_DIR}"

cat > "${SYNTH_DIR}/tiny.gtf" <<'GTF'
chr1	src	exon	100	149	.	+	.	gene_id "GENE1"; transcript_id "TX1"; gene_name "G1";
chr1	src	CDS	110	140	.	+	0	gene_id "GENE1"; transcript_id "TX1"; gene_name "G1";
chr1	src	exon	200	249	.	+	.	gene_id "GENE1"; transcript_id "TX1"; gene_name "G1";
GTF

cat > "${SYNTH_DIR}/tiny.bed" <<'BED'
chr1	120	130	tx1_read1	60	+
chr1	210	220	tx1_read2	45	-
BED

log "wrote synthetic fixtures under ${SYNTH_DIR}/"

# 3. Public GTF / BED fetches (network-dependent; skipped silently on failure).
PUBLIC_DIR="${DEST_DIR}/public"
mkdir -p "${PUBLIC_DIR}"

# Example placeholder URLs -- uncomment and adjust as needed. Kept commented
# because the maintainer community has not yet picked canonical inputs.
# fetch "https://example.org/drvizer-bench/gencode-tiny.gtf.gz" \
#       "${PUBLIC_DIR}/gencode-tiny.gtf.gz" || skip "public GTF fetch failed"

log "done. Data root: ${DEST_DIR}"
