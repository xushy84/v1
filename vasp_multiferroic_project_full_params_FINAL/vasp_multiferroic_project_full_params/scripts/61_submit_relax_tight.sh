#!/usr/bin/env bash
set -euo pipefail

RUNS_ROOT="${RUNS_ROOT:-runs}"
QSUB_BIN="${QSUB_BIN:-qsub}"
MAX_SUBMIT="${MAX_SUBMIT:-999999}"
LOG_FILE="${LOG_FILE:-submit_relax_tight.log}"

STEP_DIR="${STEP_DIR:-relax_tight}"
JOB_FILE_NAME="${JOB_FILE_NAME:-job_relax_tight.pbs}"
JOBID_FILE_NAME="pbs.jobid"
OUTCAR_COMPLETE_MARKER="General timing and accounting"

print_help() {
cat << EOF
Usage: $(basename "$0") [OPTIONS]
  -h, --help          help
  -v, --verbose       verbose
  -d, --dry-run       dry run
  -l, --log FILE      log file
  -r, --runs-root DIR runs root
  -m, --max-submit N  max submit
EOF
}

log() {
  local level="$1"; local msg="$2"
  local ts; ts=$(date +"%Y-%m-%d %H:%M:%S")
  echo "[$ts] [$level] $msg" | tee -a "$LOG_FILE"
}

is_job_running() {
  local jobid="$1"
  command -v qstat >/dev/null 2>&1 || return 1
  qstat "$jobid" >/dev/null 2>&1
}

is_job_complete() {
  local outcar="$1"
  [ -f "$outcar" ] || return 1
  grep -qi "$OUTCAR_COMPLETE_MARKER" "$outcar" 2>/dev/null
}

VERBOSE=0
DRY_RUN=0
while [[ $# -gt 0 ]]; do
  case "$1" in
    -h|--help) print_help; exit 0 ;;
    -v|--verbose) VERBOSE=1; shift ;;
    -d|--dry-run) DRY_RUN=1; shift ;;
    -l|--log) LOG_FILE="$2"; shift 2 ;;
    -r|--runs-root) RUNS_ROOT="$2"; shift 2 ;;
    -m|--max-submit) MAX_SUBMIT="$2"; shift 2 ;;
    *) log "ERROR" "Unknown arg: $1"; print_help; exit 1 ;;
  esac
done

[ -d "$RUNS_ROOT" ] || { log "ERROR" "Missing RUNS_ROOT: $RUNS_ROOT"; exit 1; }
command -v "$QSUB_BIN" >/dev/null 2>&1 || { log "ERROR" "Missing qsub: $QSUB_BIN"; exit 1; }
[[ "$MAX_SUBMIT" =~ ^[0-9]+$ ]] || { log "ERROR" "MAX_SUBMIT must be integer: $MAX_SUBMIT"; exit 1; }

total_dirs=0; complete_jobs=0; running_jobs=0; submitted_jobs=0; skipped_jobs=0
shopt -s nullglob
for step_dir in "${RUNS_ROOT}"/*/"${STEP_DIR}"; do
  total_dirs=$((total_dirs+1))
  [ -d "$step_dir" ] || { skipped_jobs=$((skipped_jobs+1)); continue; }

  step_dir=$(realpath "$step_dir")
  job_file="${step_dir}/${JOB_FILE_NAME}"
  outcar_file="${step_dir}/OUTCAR"
  jobid_file="${step_dir}/${JOBID_FILE_NAME}"

  [ -f "$job_file" ] || { skipped_jobs=$((skipped_jobs+1)); continue; }

  if is_job_complete "$outcar_file"; then
    complete_jobs=$((complete_jobs+1))
    [ $VERBOSE -eq 1 ] && log "INFO" "Done, skip: $step_dir"
    continue
  fi

  if [ -f "$jobid_file" ]; then
    jobid=$(tr -d ' \n\r' < "$jobid_file" 2>/dev/null || true)
    if [ -n "$jobid" ] && is_job_running "$jobid"; then
      running_jobs=$((running_jobs+1))
      [ $VERBOSE -eq 1 ] && log "INFO" "Running, skip: $step_dir (JobID: $jobid)"
      continue
    fi
  fi

  if [ "$submitted_jobs" -ge "$MAX_SUBMIT" ]; then
    log "INFO" "Reached MAX_SUBMIT=$MAX_SUBMIT, stop"
    break
  fi

  log "INFO" "Submit: $step_dir"
  if [ $DRY_RUN -eq 1 ]; then
    log "INFO" "[DRY] cd $step_dir && $QSUB_BIN $JOB_FILE_NAME"
    submitted_jobs=$((submitted_jobs+1))
    continue
  fi

  (
    cd "$step_dir" || exit 1
    jid=$("$QSUB_BIN" "$JOB_FILE_NAME" 2>/dev/null | awk '{print $NF}' | tr -d ' \n\r')
    [ -n "$jid" ] || { log "ERROR" "qsub no output: $step_dir"; exit 1; }
    echo "$jid" > "$JOBID_FILE_NAME"
    log "SUCCESS" "$step_dir -> $jid"
  ) || { skipped_jobs=$((skipped_jobs+1)); continue; }

  submitted_jobs=$((submitted_jobs+1))
done

log "INFO" "==== STATS ===="
log "INFO" "total dirs:   $total_dirs"
log "INFO" "completed:    $complete_jobs"
log "INFO" "running:      $running_jobs"
log "INFO" "submitted:    $submitted_jobs"
log "INFO" "skipped:      $skipped_jobs"
