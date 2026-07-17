from __future__ import annotations

import csv
import importlib.util
import json
from pathlib import Path
import shlex
import sys
import tempfile
import unittest


BENCH_ROOT = Path(__file__).resolve().parents[1]
SPEC = importlib.util.spec_from_file_location(
    "skype_bench", BENCH_ROOT / "skype_bench.py"
)
assert SPEC is not None and SPEC.loader is not None
skype_bench = importlib.util.module_from_spec(SPEC)
sys.modules[SPEC.name] = skype_bench
SPEC.loader.exec_module(skype_bench)


class SkypeBenchLimitCombinationTests(unittest.TestCase):
    def setUp(self) -> None:
        self.sample = skype_bench.Sample(
            "CELL",
            Path("/reads.fastq.gz"),
            {"hifi_nanomonsv": Path("/calls.vcf")},
        )

    def test_vcf_selection_adds_native_skype_first(self) -> None:
        cases = skype_bench.selected_cases(
            [self.sample], None, {"hifi_nanomonsv"}
        )
        self.assertEqual([method for _, method in cases], [
            "skype",
            "hifi_nanomonsv",
        ])

    def test_vcf_command_forwards_native_limit_file_to_stage_02(self) -> None:
        results_dir = Path("/benchmark-results")
        command, _, _, _, _ = skype_bench.command_for_case(
            self.sample, "hifi_nanomonsv", results_dir, 8, 2
        )
        option = next(
            argument
            for argument in command
            if argument.startswith("--option_02=")
        )
        forwarded = shlex.split(option.split("=", 1)[1])
        self.assertEqual(
            forwarded,
            [
                "--limit_combinations",
                str(
                    skype_bench.native_limit_combinations_path(
                        results_dir, self.sample.cell_line
                    )
                ),
            ],
        )

    def test_completed_vcf_must_match_native_limit_combinations(self) -> None:
        metrics = {
            "nclose_count": 1,
            "indel_count": 2,
            "denoised_relative_error": 0.5,
        }
        with tempfile.TemporaryDirectory() as temporary:
            results_dir = Path(temporary)
            native_limits_path = skype_bench.native_limit_combinations_path(
                results_dir, self.sample.cell_line
            )
            native_limits_path.parent.mkdir(parents=True)
            native_limits_path.write_text(
                json.dumps({"limit_combinations": [2, 1]}),
                encoding="utf-8",
            )
            native_record = {
                "status": "complete",
                "metrics": metrics,
                "limit_combinations": [2, 1],
            }
            vcf_record = {
                "status": "complete",
                "metrics": metrics,
                "limit_combinations": [2, 1],
            }
            status = {
                "runs": {
                    "CELL|skype": native_record,
                    "CELL|hifi_nanomonsv": vcf_record,
                }
            }

            self.assertIsNotNone(
                skype_bench.completed_case_metrics(
                    vcf_record,
                    self.sample,
                    "hifi_nanomonsv",
                    status,
                    results_dir,
                )
            )
            vcf_record["limit_combinations"] = [1, 1]
            self.assertIsNone(
                skype_bench.completed_case_metrics(
                    vcf_record,
                    self.sample,
                    "hifi_nanomonsv",
                    status,
                    results_dir,
                )
            )


class SkypeBenchNcloseMetricTests(unittest.TestCase):
    def setUp(self) -> None:
        self.temporary = tempfile.TemporaryDirectory()
        self.addCleanup(self.temporary.cleanup)
        self.root = Path(self.temporary.name)
        self.output_dir = self.root / "output"
        self.output_dir.mkdir()
        self.log_path = self.root / "pipeline.log"
        self.log_path.write_text(
            "INFO: Denoised relative error : 0.125\n",
            encoding="utf-8",
        )
        (self.output_dir / "vcf_mode_summary.tsv").write_text(
            "metric\tvalue\nused_type4_events\t3\n",
            encoding="utf-8",
        )

    def write_nclose_report(self, report_count: int) -> Path:
        report_path = self.output_dir / skype_bench.NCLOSE_REPORT_TSV
        with report_path.open("w", newline="", encoding="utf-8") as handle:
            writer = csv.DictWriter(
                handle,
                fieldnames=skype_bench.NCLOSE_REPORT_COLUMNS,
                delimiter="\t",
                lineterminator="\n",
            )
            writer.writeheader()
            for index in range(1, report_count + 1):
                writer.writerow({"nclose_id": f"SKYPE.nclose.{index}"})
        return report_path

    def test_collect_metrics_uses_nclose_report_rows(self) -> None:
        self.write_nclose_report(report_count=2)
        # The old source deliberately has a different count.  This verifies
        # that collect_metrics no longer reads nclose_nodes_index.txt.
        (self.output_dir / "nclose_nodes_index.txt").write_text(
            "1\n2\n3\n4\n5\n",
            encoding="utf-8",
        )

        metrics = skype_bench.collect_metrics(
            self.output_dir,
            self.log_path,
            vcf_mode=True,
        )

        self.assertEqual(metrics["nclose_count"], 2)
        self.assertEqual(metrics["indel_count"], 3)
        self.assertEqual(metrics["denoised_relative_error"], 0.125)

    def test_collect_metrics_counts_header_only_report_as_zero(self) -> None:
        self.write_nclose_report(report_count=0)

        metrics = skype_bench.collect_metrics(
            self.output_dir,
            self.log_path,
            vcf_mode=True,
        )

        self.assertEqual(metrics["nclose_count"], 0)

    def test_count_nclose_reports_rejects_unexpected_schema(self) -> None:
        report_path = self.output_dir / skype_bench.NCLOSE_REPORT_TSV
        report_path.write_text(
            "nclose_id\tstart_chr\nSKYPE.nclose.1\tchr1\n",
            encoding="utf-8",
        )

        with self.assertRaisesRegex(
            skype_bench.BenchError,
            "unexpected NClose report header",
        ):
            skype_bench.count_nclose_reports(report_path)

if __name__ == "__main__":
    unittest.main()
