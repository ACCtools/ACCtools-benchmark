from __future__ import annotations

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


if __name__ == "__main__":
    unittest.main()
