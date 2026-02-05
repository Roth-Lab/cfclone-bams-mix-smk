#!/usr/bin/env bash

DIR="."

OUT="bam_read_counts.tsv"

echo -e "file\treads" > "$OUT"

for bam in "$DIR"/*.bam; do
    if [[ -f "$bam" ]]; then
        count=$(samtools view -c "$bam")
        echo -e "$(basename "$bam")\t$count" >> "$OUT"
    fi
done

echo "Saved read counts to $OUT"

