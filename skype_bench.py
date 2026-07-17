#!/usr/bin/env python3
"""Benchmark native SKYPE karyotyping against the 11 caller VCF inputs.

For every row in cancer_raw_reads.csv this script runs:

* ``ACCtools-pipeline/SKYPE.py run_hifi`` directly once for the native,
  assembly-driven SKYPE result using the HS1 reference.  This run is always
  the first case and supplies its successful graph ``limit_combinations``.
* ``ACCtools-pipeline/SKYPE.py run_hifi`` directly for each of the four HiFi,
  four ONT, and three Illumina caller VCFs.  Every VCF case uses exactly the
  native graph-limit combination and fails without trying a fallback.

The requested summary is written to ``skype_bench_results/skype_bench.csv`` by
default.  Runs are logged separately and checkpointed in
``skype_bench_results/status.json``, so an interrupted invocation can be
resumed by running the same command again.

Before the first runnable case for each cell line, the shared ``30_skype`` and
``31_skype_hg38`` pipeline output directories are removed once.  Completed
cases remain available from their snapshots under the benchmark results
directory.

Metric definitions:

* nclose_count: number of report rows emitted in nclose_report.tsv
* indel_count: native type-4 indels, or VCF-mode used_type4_events
* denoised_relative_error: ``Denoised relative error`` emitted by
  23_run_nnls.py, using the CASTLE-HiFi-calibrated chromosome/run-wise TV
  target with ``lambda = 3 * noise_sigma``
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

THREAD = 16
DEPTH = 2

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
    "nclose_report.tsv",
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
    "limit_combinations.json",
)
CSV_COLUMNS = (
    "cell_line",
    "karyotype_type",
    "nclose_count",
    "indel_count",
    "denoised_relative_error",
)
STATUS_VERSION = 1
CELL_PIPELINE_OUTPUT_DIRS = ("30_skype", "31_skype_hg38")
LIMIT_COMBINATIONS_JSON = "limit_combinations.json"
NCLOSE_REPORT_TSV = "nclose_report.tsv"
NCLOSE_REPORT_COLUMNS = (
    "nclose_id",
    "start_chr",
    "start_pos",
    "start_dir",
    "end_chr",
    "end_pos",
    "end_dir",
    "nclose_cn",
    "nclose_cn_reason",
    "nclose_filter",
    "nclose_filter_reason",
    "nclose_cluster",
    "nclose_cluster_reason",
)


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


def validate_limit_combinations(
    value: object, description: str
) -> tuple[int, int]:
    if (
        not isinstance(value, (list, tuple))
        or len(value) != 2
        or any(type(item) is not int for item in value)
    ):
        raise BenchError(f"{description} must be a pair of integers")
    chr_limit, dir_limit = value
    if chr_limit < 1 or dir_limit not in (0, 1):
        raise BenchError(f"invalid {description}: {value!r}")
    return chr_limit, dir_limit


def read_limit_combinations(path: Path) -> tuple[int, int]:
    require_nonempty_file(path, "limit_combinations output")
    try:
        data = json.loads(path.read_text(encoding="utf-8"))
    except (OSError, json.JSONDecodeError) as exc:
        raise BenchError(
            f"could not read limit_combinations from {path}: {exc}"
        ) from exc
    if not isinstance(data, dict):
        raise BenchError(f"malformed limit_combinations output: {path}")
    return validate_limit_combinations(
        data.get("limit_combinations"), f"limit_combinations in {path}"
    )


def record_limit_combinations(record: object) -> tuple[int, int] | None:
    if not isinstance(record, dict):
        return None
    try:
        return validate_limit_combinations(
            record.get("limit_combinations"), "checkpoint limit_combinations"
        )
    except BenchError:
        return None


def native_limit_combinations_path(results_dir: Path, cell_line: str) -> Path:
    return (
        results_dir
        / "artifacts"
        / cell_line
        / "skype"
        / LIMIT_COMBINATIONS_JSON
    )


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
    # A legacy checkpoint has only ``relative_error``.  Treat it as pending so
    # the pipeline is rerun and actually emits the new denoised metric; deriving
    # it from the rounded legacy log value would not be possible.
    if not all(column in metrics for column in CSV_COLUMNS[2:]):
        return None
    # Runs produced before native graph limits were recorded are not
    # comparable with the new benchmark protocol and must be rerun.
    if record_limit_combinations(record) is None:
        return None
    return metrics


def native_limit_combinations_from_status(
    status: dict[str, Any], sample: Sample, results_dir: Path
) -> tuple[int, int]:
    native_record = status["runs"].get(run_key(sample.cell_line, "skype"))
    if completed_metrics(native_record) is None:
        raise BenchError(
            f"native skype prerequisite is not complete for {sample.cell_line}"
        )

    recorded = record_limit_combinations(native_record)
    if recorded is None:
        raise BenchError(
            f"native skype checkpoint has no limit_combinations for {sample.cell_line}"
        )
    artifact_path = native_limit_combinations_path(results_dir, sample.cell_line)
    artifact_value = read_limit_combinations(artifact_path)
    if artifact_value != recorded:
        raise BenchError(
            f"native limit_combinations mismatch for {sample.cell_line}: "
            f"checkpoint={recorded}, artifact={artifact_value}"
        )
    return artifact_value


def completed_case_metrics(
    record: object,
    sample: Sample,
    method: str,
    status: dict[str, Any],
    results_dir: Path,
) -> dict[str, Any] | None:
    metrics = completed_metrics(record)
    if metrics is None:
        return None
    try:
        native_limits = native_limit_combinations_from_status(
            status, sample, results_dir
        )
    except BenchError:
        return None
    if record_limit_combinations(record) != native_limits:
        return None
    return metrics


def write_summary_csv(
    output_csv: Path,
    samples: Iterable[Sample],
    status: dict[str, Any],
    results_dir: Path,
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
                metrics = completed_case_metrics(
                    runs.get(run_key(sample.cell_line, method)),
                    sample,
                    method,
                    status,
                    results_dir,
                )
                if metrics is None:
                    continue
                writer.writerow(
                    {
                        "cell_line": sample.cell_line,
                        "karyotype_type": method,
                        "nclose_count": metrics["nclose_count"],
                        "indel_count": metrics["indel_count"],
                        "denoised_relative_error": (
                            f'{float(metrics["denoised_relative_error"]):.4f}'
                        ),
                    }
                )
                rows_written += 1
    temporary.replace(output_csv)
    return rows_written


def count_nclose_reports(path: Path) -> int:
    require_nonempty_file(path, "NClose report")
    with path.open(newline="", encoding="utf-8", errors="strict") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        if tuple(reader.fieldnames or ()) != NCLOSE_REPORT_COLUMNS:
            raise BenchError(
                f"unexpected NClose report header in {path}: {reader.fieldnames}"
            )
        return sum(1 for _row in reader)


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


DENOISED_RELATIVE_ERROR_RE = re.compile(
    r"(?:^|\s)INFO:\s*Denoised relative error\s*:\s*"
    r"([0-9]+(?:\.[0-9]+)?(?:[eE][+-]?[0-9]+)?)\s*$"
)


def read_denoised_relative_error(log_path: Path) -> float:
    require_nonempty_file(log_path, "pipeline log")
    values: list[float] = []
    with log_path.open(encoding="utf-8", errors="replace") as handle:
        for line in handle:
            match = DENOISED_RELATIVE_ERROR_RE.search(line)
            if match:
                values.append(float(match.group(1)))
    if not values:
        raise BenchError(f"Denoised relative error was not found in {log_path}")
    return values[-1]


def collect_metrics(output_dir: Path, log_path: Path, vcf_mode: bool) -> dict[str, Any]:
    nclose_count = count_nclose_reports(output_dir / NCLOSE_REPORT_TSV)
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
        "denoised_relative_error": read_denoised_relative_error(log_path),
    }


def ensure_fresh_output(output_dir: Path, vcf_mode: bool, started_ns: int) -> None:
    required = [
        output_dir / NCLOSE_REPORT_TSV,
        output_dir / "karyotype.txt",
        output_dir / LIMIT_COMBINATIONS_JSON,
    ]
    if vcf_mode:
        required.append(output_dir / "vcf_mode_summary.tsv")
    else:
        required.append(output_dir / "conjoined_type4_ins_del.pkl")
    for path in required:
        require_file(path, "pipeline output")
        if path.stat().st_size == 0:
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
    needs_native_limits = any(method != "skype" for method in selected_methods)
    for sample in samples:
        if selected_cells is not None and sample.cell_line not in selected_cells:
            continue
        for method in METHOD_ORDER:
            if method in selected_methods or (
                method == "skype" and needs_native_limits
            ):
                cases.append((sample, method))
    return cases


def limit_vcf_cases(
    cases: list[tuple[Sample, str]],
    max_vcfs: int | None,
    status: dict[str, Any] | None = None,
    results_dir: Path | None = None,
    force: bool = False,
) -> list[tuple[Sample, str]]:
    if max_vcfs is None:
        return cases

    limited: list[tuple[Sample, str]] = []
    selected_vcfs: dict[str, int] = {}
    completed_runs = status["runs"] if status is not None else {}
    for sample, method in cases:
        if method == "skype":
            limited.append((sample, method))
            continue
        if (
            not force
            and status is not None
            and results_dir is not None
            and completed_case_metrics(
                completed_runs.get(run_key(sample.cell_line, method)),
                sample,
                method,
                status,
                results_dir,
            ) is not None
        ):
            limited.append((sample, method))
        elif selected_vcfs.get(sample.cell_line, 0) < max_vcfs:
            limited.append((sample, method))
            selected_vcfs[sample.cell_line] = (
                selected_vcfs.get(sample.cell_line, 0) + 1
            )
    return limited


def command_for_case(
    sample: Sample, method: str, results_dir: Path, thread: int, depth: int
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
            str(thread),
            "-d",
            str(depth),
            str(cell_root),
            str(sample.raw_read),
        ]
        output_dir = cell_root / "30_skype"
        return command, output_dir, log_path, False, artifact_dir

    vcf_path = sample.vcfs[method]
    native_limits_path = native_limit_combinations_path(
        results_dir, sample.cell_line
    )
    option_02 = shlex.join(
        ["--limit_combinations", str(native_limits_path)]
    )
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
        f"--option_02={option_02}",
        "-t",
        str(thread),
        "-d",
        str(depth),
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
    sample: Sample,
    method: str,
    results_dir: Path,
    thread: int,
    depth: int,
    expected_limit_combinations: tuple[int, int] | None = None,
) -> dict[str, Any]:
    command, output_dir, log_path, vcf_mode, snapshot_dir = command_for_case(
        sample, method, results_dir, thread, depth
    )
    if vcf_mode:
        if expected_limit_combinations is None:
            raise BenchError(
                f"native limit_combinations were not supplied for "
                f"{sample.cell_line}/{method}"
            )
        native_limits_path = native_limit_combinations_path(
            results_dir, sample.cell_line
        )
        source_limit_combinations = read_limit_combinations(native_limits_path)
        if source_limit_combinations != expected_limit_combinations:
            raise BenchError(
                f"native limit_combinations changed before "
                f"{sample.cell_line}/{method}: expected="
                f"{expected_limit_combinations}, actual={source_limit_combinations}"
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
    actual_limit_combinations = read_limit_combinations(
        output_dir / LIMIT_COMBINATIONS_JSON
    )
    if (
        expected_limit_combinations is not None
        and actual_limit_combinations != expected_limit_combinations
    ):
        raise BenchError(
            f"pipeline used unexpected limit_combinations for "
            f"{sample.cell_line}/{method}: expected={expected_limit_combinations}, "
            f"actual={actual_limit_combinations}"
        )
    metrics = collect_metrics(output_dir, log_path, vcf_mode)
    base_record.update(
        {
            "status": "complete",
            "metrics": metrics,
            "limit_combinations": list(actual_limit_combinations),
        }
    )

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
        "--max-vcfs",
        type=int,
        metavar="N",
        help=(
            "Run at most N pending VCF cases per cell line; completed "
            "checkpoints and native skype cases do not count."
        ),
    )
    parser.add_argument(
        "-t",
        "--thread",
        type=int,
        default=THREAD,
        help=f"Number of pipeline threads (default: {THREAD}).",
    )
    parser.add_argument(
        "-d",
        "--depth",
        type=int,
        default=DEPTH,
        help=f"Breakend graph depth (default: {DEPTH}).",
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
    if args.max_vcfs is not None and args.max_vcfs < 1:
        raise BenchError("--max-vcfs must be at least 1")
    if args.thread < 1:
        raise BenchError("--thread must be at least 1")
    if args.depth < 1:
        raise BenchError("--depth must be at least 1")
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

    if args.dry_run:
        cases = limit_vcf_cases(cases, args.max_vcfs)

    print(
        f"Validated {len(samples)} cell lines, 11 VCFs per cell; "
        f"selected {len(cases)} benchmark cases.",
        flush=True,
    )
    if args.dry_run:
        for sample, method in cases:
            command, output_dir, log_path, _, _ = command_for_case(
                sample, method, results_dir, args.thread, args.depth
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
        cases = limit_vcf_cases(
            cases,
            args.max_vcfs,
            status=status,
            results_dir=results_dir,
            force=args.force,
        )
        failures = 0
        cleaned_cells: set[str] = set()

        for index, (sample, method) in enumerate(cases, start=1):
            key = run_key(sample.cell_line, method)
            previous = status["runs"].get(key)
            if (
                not args.force
                and completed_case_metrics(
                    previous, sample, method, status, results_dir
                ) is not None
            ):
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
                expected_limit_combinations = None
                if method != "skype":
                    expected_limit_combinations = (
                        native_limit_combinations_from_status(
                            status, sample, results_dir
                        )
                    )
                if sample.cell_line not in cleaned_cells:
                    remove_cell_pipeline_outputs(sample.cell_line)
                    cleaned_cells.add(sample.cell_line)
                record = run_one_case(
                    sample,
                    method,
                    results_dir,
                    args.thread,
                    args.depth,
                    expected_limit_combinations,
                )
            except KeyboardInterrupt:
                status["runs"][key] = {
                    "status": "interrupted",
                    "cell_line": sample.cell_line,
                    "karyotype_type": method,
                    "finished_at": now_iso(),
                }
                atomic_write_json(status_path, status)
                write_summary_csv(output_csv, samples, status, results_dir)
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
            rows = write_summary_csv(output_csv, samples, status, results_dir)
            if record.get("status") == "complete":
                metrics = record["metrics"]
                print(
                    f"METRIC {sample.cell_line}/{method}: "
                    f"nclose={metrics['nclose_count']} "
                    f"indel={metrics['indel_count']} "
                    "denoised_relative_error="
                    f"{metrics['denoised_relative_error']:.4f} "
                    f"(CSV rows={rows})",
                    flush=True,
                )
            if failures and args.fail_fast:
                break

        rows = write_summary_csv(output_csv, samples, status, results_dir)
        selected_incomplete = [
            run_key(sample.cell_line, method)
            for sample, method in cases
            if completed_case_metrics(
                status["runs"].get(run_key(sample.cell_line, method)),
                sample,
                method,
                status,
                results_dir,
            ) is None
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
