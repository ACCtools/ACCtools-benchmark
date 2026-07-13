#!/usr/bin/env python3
"""Benchmark native SKYPE karyotyping against the 11 caller VCF inputs.

For every row in cancer_raw_reads.csv this script runs:

* ``ACCtools-pipeline/SKYPE.py run_hifi`` directly once for the native,
  assembly-driven SKYPE result using the HS1 reference.
* ``ACCtools-pipeline/SKYPE.py run_hifi`` directly for each of the four HiFi,
  four ONT, and three Illumina caller VCFs.

The requested summary is written to ``skype_bench_results/skype_bench.csv`` by
default.  Runs are logged separately and checkpointed in
``skype_bench_results/status.json``, so an interrupted invocation can be
resumed by running the same command again.

Before the first runnable case for each cell line, the shared ``30_skype`` and
``31_skype_hg38`` pipeline output directories are removed once.  Completed
cases remain available from their snapshots under the benchmark results
directory.

Metric definitions:

* nclose_count: number of entries in nclose_nodes_index.txt
* indel_count: native type-4 indels, or VCF-mode used_type4_events
* relative_error: the base ``Relative error`` emitted by 23_run_nnls.py
"""

from __future__ import annotations

import argparse
import csv
from dataclasses import dataclass
from datetime import datetime, timezone
import fcntl
import json
import os
from pathlib import Path
import pickle
import re
import shlex
import shutil
import signal
import subprocess
import sys
import time
from typing import Any, Iterable


PROJECT_ROOT = Path(__file__).resolve().parent
WORKSPACE_ROOT = PROJECT_ROOT.parent
THREAD = 50
DEPTH = 1
MAMBA_BIN_DIR = Path("/home/hyunwoo/.mamba/bin")
SKYPE_ENV_BIN = Path("/hyunwoo/.mamba/envs/skype/bin")
SKYPE_PYTHON = SKYPE_ENV_BIN / "python"
ACCTOOLS_SKYPE = WORKSPACE_ROOT / "ACCtools-pipeline" / "SKYPE.py"
DATA_ROOT = Path("/Data/hyunwoo/00_skype_run_data")
SNAPSHOT_FILES = (
    "karyotype.txt",
    "karyotype_filter.txt",
    "karyotype_cluster.txt",
    "nclose_nodes_index.txt",
    "conjoined_type4_ins_del.pkl",
    "vcf_type4_events.pkl",
    "vcf_mode_summary.json",
    "vcf_mode_summary.tsv",
    "vcf_mode_skipped_records.tsv",
    "vcf_mode_orientation_mismatches.tsv",
    "SV_benchmark_result.vcf",
    "SV_call_result.vcf",
    "SV_call_result_filter.vcf",
    "SV_call_result_cluster.vcf",
    "SKYPE_result.bed",
    "SKYPE_result_filter.bed",
    "SKYPE_result_cluster.bed",
    "total_cov.png",
    "virtual_sky.png",
    "*.paf.ppc.paf",
    "compressed_nclose_nodes_list.txt",
    "all_nclose_nodes_list.txt",
    "report.txt",
    "pipeline_mode.pkl",
)
CSV_COLUMNS = (
    "cell_line",
    "karyotype_type",
    "nclose_count",
    "indel_count",
    "relative_error",
)
STATUS_VERSION = 1
CELL_PIPELINE_OUTPUT_DIRS = ("30_skype", "31_skype_hg38")


@dataclass(frozen=True)
class VcfMethod:
    name: str
    pattern: str


@dataclass(frozen=True)
class Sample:
    cell_line: str
    raw_read: Path
    vcfs: dict[str, Path]


VCF_METHODS = (
    VcfMethod("hifi_nanomonsv", "tumor.nanomonsv.result_HIFI_vaf.vcf"),
    VcfMethod("hifi_savana", "savana_*_HIFI_vaf.vcf"),
    VcfMethod("hifi_severus", "severus_HIFI_vaf.vcf"),
    VcfMethod("hifi_sniffles2", "sniffles2_HIFI_vaf.vcf"),
    VcfMethod("ont_nanomonsv", "tumor.nanomonsv.result_ONT_vaf.vcf"),
    VcfMethod("ont_savana", "savana_*_ONT_vaf.vcf"),
    VcfMethod("ont_severus", "severus_ont_vaf.vcf"),
    VcfMethod("ont_sniffles2", "sniffles2_ONT_vaf.vcf"),
    VcfMethod("illumina_gripss", "*.gripss.filtered_vaf.vcf"),
    VcfMethod("illumina_manta", "manta_somaticSV_vaf.vcf"),
    VcfMethod("illumina_svaba", "svaba_merged.vcf"),
)
METHOD_ORDER = ("skype", *(method.name for method in VCF_METHODS))


class BenchError(RuntimeError):
    """Raised for an invalid benchmark input or result."""


class PipelineError(BenchError):
    def __init__(self, message: str, returncode: int):
        super().__init__(message)
        self.returncode = returncode


def now_iso() -> str:
    return datetime.now(timezone.utc).isoformat()


def absolute_path(path: Path) -> Path:
    return path.expanduser().resolve()


def require_nonempty_file(path: Path, description: str) -> None:
    if not path.is_file():
        raise BenchError(f"{description} not found: {path}")
    if path.stat().st_size == 0:
        raise BenchError(f"{description} is empty: {path}")


def require_file(path: Path, description: str) -> None:
    if not path.is_file():
        raise BenchError(f"{description} not found: {path}")


def read_samples(input_csv: Path, vcf_root: Path) -> list[Sample]:
    require_nonempty_file(input_csv, "input CSV")
    samples: list[Sample] = []
    seen_cells: set[str] = set()

    with input_csv.open(newline="", encoding="utf-8-sig") as handle:
        reader = csv.DictReader(handle)
        required_columns = {"cancer_prefix", "raw_read"}
        if reader.fieldnames is None or not required_columns.issubset(reader.fieldnames):
            raise BenchError(
                f"{input_csv} must contain columns: cancer_prefix, raw_read"
            )

        for line_number, row in enumerate(reader, start=2):
            cell = (row.get("cancer_prefix") or "").strip()
            raw_value = (row.get("raw_read") or "").strip()
            if not cell or not raw_value:
                raise BenchError(
                    f"{input_csv}:{line_number}: cancer_prefix/raw_read is empty"
                )
            if cell in seen_cells:
                raise BenchError(f"duplicate cell line in {input_csv}: {cell}")
            seen_cells.add(cell)

            raw_read = absolute_path(Path(raw_value))
            require_nonempty_file(raw_read, f"raw read for {cell}")

            vcf_dir = vcf_root / cell / cell
            if not vcf_dir.is_dir():
                raise BenchError(f"VCF directory not found for {cell}: {vcf_dir}")

            vcfs: dict[str, Path] = {}
            used_paths: set[Path] = set()
            for method in VCF_METHODS:
                matches = sorted(path.resolve() for path in vcf_dir.glob(method.pattern))
                matches = [path for path in matches if path.is_file()]
                if len(matches) != 1:
                    rendered = ", ".join(str(path) for path in matches) or "none"
                    raise BenchError(
                        f"{cell}/{method.name}: expected one VCF matching "
                        f"{method.pattern!r}, found {len(matches)} ({rendered})"
                    )
                path = matches[0]
                require_nonempty_file(path, f"VCF for {cell}/{method.name}")
                if path in used_paths:
                    raise BenchError(f"VCF matched more than one method: {path}")
                used_paths.add(path)
                vcfs[method.name] = path

            if len(vcfs) != 11:
                raise BenchError(f"{cell}: expected 11 VCFs, found {len(vcfs)}")
            samples.append(Sample(cell, raw_read, vcfs))

    if not samples:
        raise BenchError(f"no samples found in {input_csv}")
    return samples


def load_status(path: Path) -> dict[str, Any]:
    if not path.exists():
        return {"version": STATUS_VERSION, "runs": {}}
    try:
        status = json.loads(path.read_text(encoding="utf-8"))
    except (OSError, json.JSONDecodeError) as exc:
        raise BenchError(f"could not read status file {path}: {exc}") from exc
    if status.get("version") != STATUS_VERSION or not isinstance(
        status.get("runs"), dict
    ):
        raise BenchError(f"unsupported or malformed status file: {path}")
    return status


def atomic_write_json(path: Path, data: dict[str, Any]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    temporary = path.with_name(f".{path.name}.tmp-{os.getpid()}")
    with temporary.open("w", encoding="utf-8") as handle:
        json.dump(data, handle, indent=2, sort_keys=True)
        handle.write("\n")
    temporary.replace(path)


def run_key(cell_line: str, method: str) -> str:
    return f"{cell_line}|{method}"


def completed_metrics(record: object) -> dict[str, Any] | None:
    if not isinstance(record, dict) or record.get("status") != "complete":
        return None
    metrics = record.get("metrics")
    if not isinstance(metrics, dict):
        return None
    if not all(column in metrics for column in CSV_COLUMNS[2:]):
        return None
    return metrics


def write_summary_csv(
    output_csv: Path, samples: Iterable[Sample], status: dict[str, Any]
) -> int:
    output_csv.parent.mkdir(parents=True, exist_ok=True)
    temporary = output_csv.with_name(f".{output_csv.name}.tmp-{os.getpid()}")
    rows_written = 0
    runs = status["runs"]

    with temporary.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=CSV_COLUMNS)
        writer.writeheader()
        for sample in samples:
            for method in METHOD_ORDER:
                metrics = completed_metrics(runs.get(run_key(sample.cell_line, method)))
                if metrics is None:
                    continue
                writer.writerow(
                    {
                        "cell_line": sample.cell_line,
                        "karyotype_type": method,
                        "nclose_count": metrics["nclose_count"],
                        "indel_count": metrics["indel_count"],
                        "relative_error": f'{float(metrics["relative_error"]):.4f}',
                    }
                )
                rows_written += 1
    temporary.replace(output_csv)
    return rows_written


def count_nonempty_lines(path: Path) -> int:
    require_file(path, "NClose index")
    with path.open(encoding="utf-8", errors="strict") as handle:
        return sum(1 for line in handle if line.strip())


def read_metric_tsv(path: Path) -> dict[str, str]:
    require_nonempty_file(path, "VCF mode summary")
    metrics: dict[str, str] = {}
    with path.open(newline="", encoding="utf-8") as handle:
        reader = csv.reader(handle, delimiter="\t")
        header = next(reader, None)
        if header != ["metric", "value"]:
            raise BenchError(f"unexpected VCF summary header in {path}: {header}")
        for row in reader:
            if len(row) == 2:
                metrics[row[0]] = row[1]
    return metrics


def count_native_indels(path: Path) -> int:
    require_nonempty_file(path, "native type-4 indel data")
    try:
        with path.open("rb") as handle:
            data = pickle.load(handle)
    except (OSError, pickle.UnpicklingError, EOFError) as exc:
        raise BenchError(f"could not read native indel data {path}: {exc}") from exc
    if not isinstance(data, (tuple, list)) or len(data) != 2:
        raise BenchError(f"unexpected native indel data structure in {path}")
    if not all(hasattr(group, "__len__") for group in data):
        raise BenchError(f"native indel groups are not countable in {path}")
    return sum(len(group) for group in data)


RELATIVE_ERROR_RE = re.compile(
    r"(?:^|\s)INFO:\s*Relative error\s*:\s*"
    r"([0-9]+(?:\.[0-9]+)?(?:[eE][+-]?[0-9]+)?)\s*$"
)


def read_relative_error(log_path: Path) -> float:
    require_nonempty_file(log_path, "pipeline log")
    values: list[float] = []
    with log_path.open(encoding="utf-8", errors="replace") as handle:
        for line in handle:
            match = RELATIVE_ERROR_RE.search(line)
            if match:
                values.append(float(match.group(1)))
    if not values:
        raise BenchError(f"base Relative error was not found in {log_path}")
    return values[-1]


def collect_metrics(output_dir: Path, log_path: Path, vcf_mode: bool) -> dict[str, Any]:
    nclose_count = count_nonempty_lines(output_dir / "nclose_nodes_index.txt")
    if vcf_mode:
        summary = read_metric_tsv(output_dir / "vcf_mode_summary.tsv")
        value = summary.get("used_type4_events")
        if value is None:
            raise BenchError(
                f"used_type4_events missing from {output_dir / 'vcf_mode_summary.tsv'}"
            )
        try:
            indel_count = int(value)
        except ValueError as exc:
            raise BenchError(f"invalid used_type4_events value: {value!r}") from exc
    else:
        indel_count = count_native_indels(
            output_dir / "conjoined_type4_ins_del.pkl"
        )
    return {
        "nclose_count": nclose_count,
        "indel_count": indel_count,
        "relative_error": read_relative_error(log_path),
    }


def ensure_fresh_output(output_dir: Path, vcf_mode: bool, started_ns: int) -> None:
    required = [output_dir / "nclose_nodes_index.txt", output_dir / "karyotype.txt"]
    if vcf_mode:
        required.append(output_dir / "vcf_mode_summary.tsv")
    else:
        required.append(output_dir / "conjoined_type4_ins_del.pkl")
    for path in required:
        require_file(path, "pipeline output")
        if path.name != "nclose_nodes_index.txt" and path.stat().st_size == 0:
            raise BenchError(f"pipeline output is empty: {path}")
        # Allow two seconds for filesystems with coarse timestamp precision.
        if path.stat().st_mtime_ns + 2_000_000_000 < started_ns:
            raise BenchError(f"pipeline output was not refreshed: {path}")


def unique_backup_path(path: Path) -> Path:
    stamp = datetime.now().strftime("%Y%m%d-%H%M%S")
    candidate = path.with_name(f"{path.name}.previous-{stamp}")
    counter = 1
    while candidate.exists():
        candidate = path.with_name(f"{path.name}.previous-{stamp}-{counter}")
        counter += 1
    return candidate


def snapshot_output(source: Path, destination: Path, metadata: dict[str, Any]) -> None:
    if destination.exists():
        destination.rename(unique_backup_path(destination))
    destination.mkdir(parents=True, exist_ok=False)
    copied: list[str] = []
    for pattern in SNAPSHOT_FILES:
        for source_path in sorted(source.glob(pattern)):
            if source_path.is_file():
                shutil.copy2(source_path, destination / source_path.name)
                copied.append(source_path.name)
    if not copied:
        raise BenchError(f"no benchmark artifacts found to snapshot in {source}")
    manifest = dict(metadata)
    manifest["source_output_dir"] = str(source)
    manifest["copied_files"] = copied
    atomic_write_json(destination / "manifest.json", manifest)


def tail_log(path: Path, lines: int = 20) -> str:
    try:
        content = path.read_text(encoding="utf-8", errors="replace").splitlines()
    except OSError:
        return ""
    return "\n".join(content[-lines:])


def pipeline_environment() -> dict[str, str]:
    environment = os.environ.copy()
    current_path = environment.get("PATH", "")
    environment["PATH"] = f"{SKYPE_ENV_BIN}:{MAMBA_BIN_DIR}:{current_path}"
    environment.setdefault("MAMBA_ROOT_PREFIX", "/home/hyunwoo/.mamba")
    return environment


def run_pipeline(command: list[str], log_path: Path) -> tuple[int, int, str, str]:
    log_path.parent.mkdir(parents=True, exist_ok=True)
    started_at = now_iso()
    started_ns = time.time_ns()
    start_monotonic = time.monotonic()
    print(f"START {started_at}  {shlex.join(command)}", flush=True)

    with log_path.open("w", encoding="utf-8") as log_handle:
        log_handle.write(f"# started_at={started_at}\n")
        log_handle.write(f"# cwd={WORKSPACE_ROOT}\n")
        log_handle.write(f"# command={shlex.join(command)}\n")
        log_handle.flush()
        process = subprocess.Popen(
            command,
            cwd=WORKSPACE_ROOT,
            env=pipeline_environment(),
            stdout=log_handle,
            stderr=subprocess.STDOUT,
            start_new_session=True,
        )
        try:
            returncode = process.wait()
        except KeyboardInterrupt:
            os.killpg(process.pid, signal.SIGTERM)
            try:
                process.wait(timeout=30)
            except subprocess.TimeoutExpired:
                os.killpg(process.pid, signal.SIGKILL)
                process.wait()
            raise

    finished_at = now_iso()
    elapsed = time.monotonic() - start_monotonic
    print(
        f"DONE  {finished_at}  exit={returncode} elapsed={elapsed / 60:.1f}m",
        flush=True,
    )
    return returncode, started_ns, started_at, finished_at


def selected_cases(
    samples: list[Sample], selected_cells: set[str] | None, selected_methods: set[str]
) -> list[tuple[Sample, str]]:
    cases: list[tuple[Sample, str]] = []
    for sample in samples:
        if selected_cells is not None and sample.cell_line not in selected_cells:
            continue
        for method in METHOD_ORDER:
            if method in selected_methods:
                cases.append((sample, method))
    return cases


def command_for_case(
    sample: Sample, method: str, results_dir: Path
) -> tuple[list[str], Path, Path, bool, Path | None]:
    log_path = results_dir / "logs" / sample.cell_line / f"{method}.log"
    artifact_dir = results_dir / "artifacts" / sample.cell_line / method
    if method == "skype":
        cell_root = DATA_ROOT / sample.cell_line
        command = [
            str(SKYPE_PYTHON),
            str(ACCTOOLS_SKYPE),
            "run_hifi",
            "--dependency_loc",
            "deps",
            "-p",
            sample.cell_line,
            "--reference",
            "hs1",
            "--option_02=--variant_mode",
            "--skype_force",
            "-t",
            str(THREAD),
            "-d",
            str(DEPTH),
            str(cell_root),
            str(sample.raw_read),
        ]
        output_dir = cell_root / "30_skype"
        return command, output_dir, log_path, False, artifact_dir

    vcf_path = sample.vcfs[method]
    command = [
        str(SKYPE_PYTHON),
        str(ACCTOOLS_SKYPE),
        "run_hifi",
        "--dependency_loc",
        "deps",
        "-p",
        sample.cell_line,
        "--reference",
        "hg38",
        "--benchmark_vcf_loc",
        str(vcf_path),
        "-t",
        str(THREAD),
        "-d",
        str(DEPTH),
        "--skype_force",
        str(DATA_ROOT / sample.cell_line),
        str(sample.raw_read),
    ]
    output_dir = DATA_ROOT / sample.cell_line / "31_skype_hg38"
    return command, output_dir, log_path, True, artifact_dir


def remove_cell_pipeline_outputs(cell_line: str) -> None:
    cell_root = DATA_ROOT / cell_line
    removed: list[Path] = []
    for directory_name in CELL_PIPELINE_OUTPUT_DIRS:
        output_dir = cell_root / directory_name
        if output_dir.is_symlink():
            raise BenchError(f"refusing to remove symlinked pipeline output: {output_dir}")
        if not output_dir.exists():
            continue
        if not output_dir.is_dir():
            raise BenchError(f"pipeline output is not a directory: {output_dir}")
        shutil.rmtree(output_dir)
        removed.append(output_dir)

    rendered = ", ".join(str(path) for path in removed) or "nothing to remove"
    print(f"CLEAN {cell_line}: {rendered}", flush=True)


def run_one_case(
    sample: Sample, method: str, results_dir: Path
) -> dict[str, Any]:
    command, output_dir, log_path, vcf_mode, snapshot_dir = command_for_case(
        sample, method, results_dir
    )
    returncode, started_ns, started_at, finished_at = run_pipeline(command, log_path)
    base_record: dict[str, Any] = {
        "cell_line": sample.cell_line,
        "karyotype_type": method,
        "command": command,
        "log": str(log_path),
        "started_at": started_at,
        "finished_at": finished_at,
        "returncode": returncode,
    }
    if method != "skype":
        base_record["vcf"] = str(sample.vcfs[method])

    if returncode != 0:
        detail = tail_log(log_path)
        if detail:
            base_record["log_tail"] = detail
        raise PipelineError(
            f"pipeline failed for {sample.cell_line}/{method}; see {log_path}",
            returncode,
        )

    ensure_fresh_output(output_dir, vcf_mode, started_ns)
    metrics = collect_metrics(output_dir, log_path, vcf_mode)
    base_record.update({"status": "complete", "metrics": metrics})

    if snapshot_dir is not None:
        snapshot_output(output_dir, snapshot_dir, base_record)
        base_record["artifacts"] = str(snapshot_dir)
    else:
        base_record["artifacts"] = str(output_dir)
    return base_record


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Benchmark native SKYPE against 11 HiFi/ONT/Illumina VCF callers."
    )
    parser.add_argument(
        "--input",
        type=Path,
        default=PROJECT_ROOT / "cancer_raw_reads.csv",
        help="Input CSV with cancer_prefix and raw_read columns.",
    )
    parser.add_argument(
        "--vcf-root",
        type=Path,
        default=WORKSPACE_ROOT / "variant_calls_and_benchmarks",
        help="Root containing CELL_LINE/CELL_LINE/*.vcf.",
    )
    parser.add_argument(
        "--output",
        type=Path,
        help=(
            "Summary CSV path. Defaults to "
            "<results-dir>/skype_bench.csv."
        ),
    )
    parser.add_argument(
        "--results-dir",
        type=Path,
        default=WORKSPACE_ROOT / "skype_bench_results",
        help="Logs, checkpoints, and compact VCF artifacts.",
    )
    parser.add_argument(
        "--cell",
        action="append",
        dest="cells",
        help="Run only this cell line; repeat to select multiple cells.",
    )
    parser.add_argument(
        "--method",
        action="append",
        choices=METHOD_ORDER,
        help="Run only this karyotype type; repeat to select multiple types.",
    )
    parser.add_argument(
        "--force",
        action="store_true",
        help="Rerun completed selected cases; previous artifacts are preserved.",
    )
    parser.add_argument(
        "--fail-fast",
        action="store_true",
        help="Stop at the first failed pipeline instead of continuing.",
    )
    parser.add_argument(
        "--dry-run",
        action="store_true",
        help="Validate all inputs and print commands without running pipelines.",
    )
    return parser.parse_args()


def main() -> int:
    args = parse_args()
    input_csv = absolute_path(args.input)
    vcf_root = absolute_path(args.vcf_root)
    results_dir = absolute_path(args.results_dir)
    output_csv = (
        absolute_path(args.output)
        if args.output is not None
        else results_dir / "skype_bench.csv"
    )

    require_nonempty_file(ACCTOOLS_SKYPE, "ACCtools SKYPE runner")
    require_nonempty_file(SKYPE_PYTHON, "skype environment Python")
    require_nonempty_file(MAMBA_BIN_DIR / "mamba", "mamba executable")

    samples = read_samples(input_csv, vcf_root)
    known_cells = {sample.cell_line for sample in samples}
    selected_cells = set(args.cells) if args.cells else None
    if selected_cells is not None:
        unknown = sorted(selected_cells - known_cells)
        if unknown:
            raise BenchError(f"--cell not present in {input_csv}: {', '.join(unknown)}")
    selected_methods = set(args.method or METHOD_ORDER)
    cases = selected_cases(samples, selected_cells, selected_methods)
    if not cases:
        raise BenchError("no benchmark cases selected")

    print(
        f"Validated {len(samples)} cell lines, 11 VCFs per cell; "
        f"selected {len(cases)} benchmark cases.",
        flush=True,
    )
    if args.dry_run:
        for sample, method in cases:
            command, output_dir, log_path, _, _ = command_for_case(
                sample, method, results_dir
            )
            print(f"{sample.cell_line}\t{method}\t{shlex.join(command)}")
            print(f"  output={output_dir}")
            print(f"  log={log_path}")
        return 0

    results_dir.mkdir(parents=True, exist_ok=True)
    lock_path = results_dir / ".lock"
    with lock_path.open("a+", encoding="utf-8") as lock_handle:
        try:
            fcntl.flock(lock_handle, fcntl.LOCK_EX | fcntl.LOCK_NB)
        except BlockingIOError as exc:
            raise BenchError(
                f"another skype_bench.py process holds the lock: {lock_path}"
            ) from exc

        status_path = results_dir / "status.json"
        status = load_status(status_path)
        failures = 0
        cleaned_cells: set[str] = set()

        for index, (sample, method) in enumerate(cases, start=1):
            key = run_key(sample.cell_line, method)
            previous = status["runs"].get(key)
            if not args.force and completed_metrics(previous) is not None:
                print(
                    f"SKIP  [{index}/{len(cases)}] {sample.cell_line}/{method} "
                    "(already complete)",
                    flush=True,
                )
                continue

            print(
                f"CASE  [{index}/{len(cases)}] {sample.cell_line}/{method}",
                flush=True,
            )
            try:
                if sample.cell_line not in cleaned_cells:
                    remove_cell_pipeline_outputs(sample.cell_line)
                    cleaned_cells.add(sample.cell_line)
                record = run_one_case(sample, method, results_dir)
            except KeyboardInterrupt:
                status["runs"][key] = {
                    "status": "interrupted",
                    "cell_line": sample.cell_line,
                    "karyotype_type": method,
                    "finished_at": now_iso(),
                }
                atomic_write_json(status_path, status)
                write_summary_csv(output_csv, samples, status)
                raise
            except Exception as exc:
                failures += 1
                returncode = exc.returncode if isinstance(exc, PipelineError) else None
                record = {
                    "status": "failed",
                    "cell_line": sample.cell_line,
                    "karyotype_type": method,
                    "finished_at": now_iso(),
                    "error": str(exc),
                }
                if returncode is not None:
                    record["returncode"] = returncode
                print(f"FAIL  {sample.cell_line}/{method}: {exc}", file=sys.stderr)
            status["runs"][key] = record
            atomic_write_json(status_path, status)
            rows = write_summary_csv(output_csv, samples, status)
            if record.get("status") == "complete":
                metrics = record["metrics"]
                print(
                    f"METRIC {sample.cell_line}/{method}: "
                    f"nclose={metrics['nclose_count']} "
                    f"indel={metrics['indel_count']} "
                    f"relative_error={metrics['relative_error']:.4f} "
                    f"(CSV rows={rows})",
                    flush=True,
                )
            if failures and args.fail_fast:
                break

        rows = write_summary_csv(output_csv, samples, status)
        selected_incomplete = [
            run_key(sample.cell_line, method)
            for sample, method in cases
            if completed_metrics(status["runs"].get(run_key(sample.cell_line, method)))
            is None
        ]
        print(f"Summary: {output_csv} ({rows} completed rows)", flush=True)
        print(f"Status:  {status_path}", flush=True)
        if selected_incomplete:
            print(
                "Incomplete selected cases: " + ", ".join(selected_incomplete),
                file=sys.stderr,
            )
            return 1
    return 0


if __name__ == "__main__":
    try:
        raise SystemExit(main())
    except BenchError as exc:
        print(f"error: {exc}", file=sys.stderr)
        raise SystemExit(1) from exc
