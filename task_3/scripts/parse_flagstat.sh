#!/usr/bin/env bash

FLAGSTAT_FILE="$1"

if [ -z "$FLAGSTAT_FILE" ]; then
  echo "Usage: $0 <flagstat_file>"
  exit 1
fi

MAPPED_PERCENT=$(grep " mapped (" "$FLAGSTAT_FILE" | head -1 | sed -E 's/.*\(([0-9.]+)% :.*/\1/')

echo "Mapped reads: ${MAPPED_PERCENT}%"

awk -v p="$MAPPED_PERCENT" 'BEGIN {
  if (p > 90) {
    print "OK"
  } else {
    print "not OK"
  }
}'