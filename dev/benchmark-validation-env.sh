#!/usr/bin/env bash
# TEMPORARY validation-only (#164); removed before merge.
# On push runs, apply allowlisted OPAL_BENCHMARK_* lines from the pushed commit
# message (VALIDATION_MESSAGE), then run the measurement command unchanged.
set -euo pipefail
job=$1
shift
if [ "${GITHUB_EVENT_NAME:-}" = push ]; then
  case $job in
    memory) allowed='UPDATE|INJECT_MEMORY|INJECT_CASES' ;;
    timing-ab) allowed='TIMING_PROFILE|ROUNDS|INJECT_SLOWDOWN|INJECT_FIT_START|INJECT_CASES|INJECT_VERSION' ;;
    instructions) allowed='UPDATE|INJECT_SLOWDOWN|INJECT_CASES' ;;
    *) echo "Unknown benchmark job: $job" >&2; exit 1 ;;
  esac
  label=$(printf '%s\n' "${VALIDATION_MESSAGE:-}" | head -n 1 | tr -cd 'A-Za-z0-9 ._:,()[]-')
  export OPAL_BENCHMARK_VALIDATION_LABEL="$label"
  echo "Validation run: $label" >> "$GITHUB_STEP_SUMMARY"
  while IFS= read -r line; do
    line=${line%$'\r'}
    if [[ $line =~ ^OPAL_BENCHMARK_($allowed)=([A-Za-z0-9._,-]*)$ ]]; then
      export "OPAL_BENCHMARK_${BASH_REMATCH[1]}=${BASH_REMATCH[2]}"
      echo "Validation setting: \`$line\`" >> "$GITHUB_STEP_SUMMARY"
    fi
  done <<< "${VALIDATION_MESSAGE:-}"
fi
exec "$@"
