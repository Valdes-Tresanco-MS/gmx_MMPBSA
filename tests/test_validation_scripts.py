import os
import json
import subprocess
import sys
import tempfile
import unittest
from pathlib import Path

from scripts.validation.common import (
    Case,
    _generated_example_files,
    apply_legacy_settings,
    compare_csv,
    csv_comparison_for_cases,
    environment_for,
    load_cases,
    load_comparison_policies,
    resolve_selectors,
)
from scripts.validation.run_api_validation import discover


class ValidationHarnessTest(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.manifest, cls.cases = load_cases()

    def test_all_selector_resolves_unique_manifest_cases(self):
        selected = resolve_selectors(self.manifest, ["101", "11", "gbnsr6"])
        self.assertEqual(selected, self.manifest["suites"]["all"]["tests"])

    def test_generated_file_filter_preserves_tracked_api_fixture(self):
        ignored = _generated_example_files(
            "examples",
            [
                "FINAL_RESULTS_MMPBSA.dat",
                "gmx_MMPBSA.log",
                "_GMXMMPBSA_complex.pdb",
                "_GMXMMPBSA_COM_FIXED.pdb",
                "mmpbsa.in",
            ],
        )
        self.assertEqual(
            ignored,
            {"FINAL_RESULTS_MMPBSA.dat", "gmx_MMPBSA.log", "_GMXMMPBSA_complex.pdb"},
        )

    def test_legacy_settings_patch_only_changes_copied_input(self):
        case = Case(
            id=1,
            name="synthetic",
            path="synthetic",
            workdir="synthetic",
            input="mmpbsa.in",
            executable="gmx_MMPBSA",
            command_args=(),
            expected_outputs=(),
            slow=False,
            requires=(),
        )
        with tempfile.TemporaryDirectory() as tmpdir:
            root = Path(tmpdir)
            workdir = root / "synthetic"
            workdir.mkdir()
            input_path = workdir / "mmpbsa.in"
            input_path.write_text(
                "&general\nPBRadii=4,\n/\n&gb\nigb=8,\n/\n&pb\nexdi=78.5,\n/\n",
                encoding="utf-8",
            )
            apply_legacy_settings(root, case)
            text = input_path.read_text(encoding="utf-8")
            self.assertIn("PBRadii=3", text)
            self.assertIn("igb=5", text)
            self.assertIn("exdi=80", text)

    def test_environment_for_exposes_selected_ambertools_data_roots(self):
        python = Path("/opt/conda/envs/test/bin/python")
        with tempfile.TemporaryDirectory() as tmpdir:
            root = Path(tmpdir)
            env = environment_for(python, root)
        self.assertEqual(env["AMBERHOME"], "/opt/conda/envs/test")
        self.assertEqual(env["DFTB_PREFIX"], "/opt/conda/envs/test/dat/slko")

    def test_csv_comparison_uses_numeric_tolerance(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            root = Path(tmpdir)
            left = root / "left.csv"
            right = root / "right.csv"
            left.write_text("Frame #,TOTAL\n1,-10.0000\n", encoding="utf-8")
            right.write_text("Frame #,TOTAL\n1,-10.0005\n", encoding="utf-8")
            self.assertEqual(compare_csv(left, right, atol=1e-3, rtol=1e-6), [])
            self.assertTrue(compare_csv(left, right, atol=1e-5, rtol=1e-6))

    def test_csv_comparison_can_scope_to_delta_energy_terms(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            root = Path(tmpdir)
            left = root / "left.csv"
            right = root / "right.csv"
            left.write_text(
                "Complex Energy Terms\nFrame #,CMAP,TOTAL\n1,0.0,-10.0\n\n"
                "Delta Energy Terms\nFrame #,CMAP,TOTAL\n1,0.0,-2.0\n",
                encoding="utf-8",
            )
            right.write_text(
                "Complex Energy Terms\nFrame #,CMAP,TOTAL\n1,-40.0,-50.0\n\n"
                "Delta Energy Terms\nFrame #,CMAP,TOTAL\n1,0.0,-2.0\n",
                encoding="utf-8",
            )
            self.assertTrue(compare_csv(left, right, atol=1e-3, rtol=1e-6))
            self.assertEqual(
                compare_csv(
                    left,
                    right,
                    atol=1e-3,
                    rtol=1e-6,
                    section="Delta Energy Terms",
                ),
                [],
            )

    def test_comparison_policies_classify_scoped_and_unsupported_cases(self):
        policies = load_comparison_policies()
        self.assertEqual(policies[10]["csv_scope"], "delta")
        self.assertEqual(policies[25]["baseline"]["status"], "unsupported")

        with tempfile.TemporaryDirectory() as tmpdir:
            root = Path(tmpdir)
            current_root = root / "current"
            baseline_root = root / "baseline"
            for run_root in (current_root, baseline_root):
                (run_root / "cases").mkdir(parents=True)
                (run_root / "examples" / "case").mkdir(parents=True)
            current_log = current_root / "cases" / "10.log"
            baseline_log = baseline_root / "cases" / "10.log"
            current_log.write_text("ok", encoding="utf-8")
            baseline_log.write_text("ok", encoding="utf-8")
            csv = (
                "Complex Energy Terms\nFrame #,CMAP,TOTAL\n1,0.0,-10.0\n\n"
                "Delta Energy Terms\nFrame #,CMAP,TOTAL\n1,0.0,-2.0\n"
            )
            baseline_csv = csv.replace("1,0.0,-10.0", "1,-40.0,-50.0")
            (current_root / "examples" / "case" / "FINAL_RESULTS_MMPBSA.csv").write_text(csv)
            (baseline_root / "examples" / "case" / "FINAL_RESULTS_MMPBSA.csv").write_text(baseline_csv)
            record = {
                "id": 10,
                "status": "PASS",
                "log": str(current_log),
                "workdir": "case",
                "expected_outputs": ["FINAL_RESULTS_MMPBSA.csv"],
            }
            baseline_record = dict(record, log=str(baseline_log))
            comparison = csv_comparison_for_cases(
                {"cases": [record]},
                {"cases": [baseline_record]},
                atol=1e-3,
                rtol=1e-6,
                policies=policies,
            )["comparisons"][0]
            self.assertEqual(comparison["status"], "PASS")
            self.assertEqual(comparison["raw_status"], "EXPECTED-DIFFERENCE")
            self.assertEqual(comparison["comparison_scope"], "Delta Energy Terms")

            unsupported_current_log = current_root / "cases" / "25.log"
            unsupported_baseline_log = baseline_root / "cases" / "25.log"
            unsupported_current_log.write_text("ok", encoding="utf-8")
            unsupported_baseline_log.write_text(
                "ImportError: cannot import name 'gmxmmpbsa_amber'",
                encoding="utf-8",
            )
            unsupported = csv_comparison_for_cases(
                {
                    "cases": [
                        dict(record, id=25, log=str(unsupported_current_log)),
                    ]
                },
                {
                    "cases": [
                        dict(
                            baseline_record,
                            id=25,
                            log=str(unsupported_baseline_log),
                            status="FAIL",
                        ),
                    ]
                },
                atol=1e-3,
                rtol=1e-6,
                policies=policies,
            )["comparisons"][0]
            self.assertEqual(unsupported["status"], "BASELINE-UNSUPPORTED")

    def test_api_validation_checks_a_compact_result_in_an_external_copy(self):
        fixture = Path(__file__).resolve().parents[1] / "examples" / "API" / "COMPACT_MMXSA_RESULTS.mmxsa"
        # A legacy test leaves the process in a temporary directory that has
        # already been removed.  API loading deliberately resolves relative
        # paths while entering the result directory, so restore a live cwd
        # before exercising the external-copy harness.
        os.chdir(fixture.parents[2])
        with tempfile.TemporaryDirectory() as tmpdir:
            result = Path(tmpdir) / fixture.name
            result.write_bytes(fixture.read_bytes())
            self.assertEqual(discover([result]), [result.resolve()])
            report = Path(tmpdir) / "api-report.json"
            subprocess.run(
                [sys.executable, str(Path(__file__).resolve().parents[1] / "scripts" / "validation" / "run_api_validation.py"),
                 str(result), "--output", str(report)],
                cwd=fixture.parents[2], check=True, capture_output=True, text=True,
            )
            record = json.loads(report.read_text(encoding="utf-8"))["results"][0]
            self.assertEqual(record["status"], "PASS")
            self.assertTrue(record["checks"]["get_energy"])
            self.assertTrue(record["checks"]["get_ana_data_inmemory"])


if __name__ == "__main__":
    unittest.main()
