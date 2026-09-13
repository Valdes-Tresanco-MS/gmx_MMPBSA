import unittest
from pathlib import Path
from tempfile import TemporaryDirectory

from GMXMMPBSA.test_manifest import (
    ManifestError,
    build_help_text,
    load_manifest,
    verify_outputs,
)


class TestManifest(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.manifest = load_manifest()

    def test_load_manifest_has_all_tests(self):
        self.assertEqual(set(self.manifest.tests), set(range(3, 27)) - {11})

    def test_suite_membership(self):
        self.assertEqual(self.manifest.suites['fast']['tests'], [3, 4, 5, 7, 12, 13, 14, 15])
        self.assertEqual(self.manifest.suites['minimal']['tests'],
                         [3, 4, 5, 6, 7, 8, 10, 12, 13, 14, 15])
        self.assertEqual(
            self.manifest.suites['all']['tests'],
            [test_id for test_id in range(3, 27) if test_id != 11],
        )

    def test_resolve_suite_and_alias(self):
        all_tests = [test_id for test_id in range(3, 27) if test_id != 11]
        self.assertEqual(self.manifest.resolve_test_ids(['0']), all_tests)
        self.assertEqual(self.manifest.resolve_test_ids(['101']), all_tests)
        self.assertEqual(self.manifest.resolve_test_ids(['2']), [3, 4, 5, 7, 12, 13, 14, 15])
        self.assertEqual(self.manifest.resolve_test_ids(['11']), [6])
        self.assertEqual(self.manifest.resolve_test_ids(['6', '11']), [6])
        self.assertEqual(self.manifest.resolve_test_ids(['gbnsr6']), [24])
        self.assertEqual(self.manifest.resolve_test_ids(['3', 'gbnsr6']), [3, 24])

    def test_invalid_selector(self):
        with self.assertRaises(ManifestError):
            self.manifest.resolve_test_ids(['not-a-test'])

    def test_build_help_text_contains_latest_test(self):
        help_text = build_help_text()
        self.assertIn('* 26', help_text)
        self.assertIn('* 0      23', help_text)
        self.assertIn('Legacy alias for test 6', help_text)
        self.assertIn('* 8    x | 10  Metalloprotein-ligand', help_text)
        self.assertIn('* 9    x | 10  Multicomponent system (Comp_receptor)', help_text)
        self.assertIn('* 15   . | 10  Interaction Entropy approximation', help_text)
        self.assertIn('* 17   x | 10  Entropy calculation using Normal Mode approximation', help_text)

    def test_legacy_membrane_selector_remains_a_valid_choice(self):
        self.assertEqual(self.manifest.all_valid_choices().count('11'), 1)

    def test_all_tests_have_command_args_and_outputs(self):
        for test_id, test in self.manifest.tests.items():
            self.assertTrue(test.command_args, msg=f'test {test_id} missing command_args')
            self.assertTrue(test.expected_outputs, msg=f'test {test_id} missing expected_outputs')
            self.assertIn(test.executable, {'gmx_MMPBSA', 'amber_MMPBSA'})

    def test_decomposition_outputs(self):
        test = self.manifest.get_test(14)
        self.assertIn('FINAL_DECOMP_MMPBSA.dat', test.expected_outputs)
        self.assertIn('FINAL_DECOMP_MMPBSA.csv', test.expected_outputs)

    def test_verify_outputs(self):
        with TemporaryDirectory() as tmpdir:
            work_dir = Path(tmpdir)
            (work_dir / 'FINAL_RESULTS_MMPBSA.dat').write_text('ok')
            missing = verify_outputs(work_dir, ['FINAL_RESULTS_MMPBSA.dat', 'FINAL_RESULTS_MMPBSA.csv'])
            self.assertEqual(missing, ['FINAL_RESULTS_MMPBSA.csv'])


def resolve_test_ids(selectors):
    return load_manifest().resolve_test_ids(selectors)


class TestManifestCommandArgs(unittest.TestCase):
    def test_explicit_waters_uses_self_contained_paths(self):
        test = load_manifest().get_test(26)
        self.assertNotIn('../', ' '.join(test.command_args))

        example_dir = Path(__file__).resolve().parents[1] / 'examples' / test.workdir
        for option in ('-i', '-cs', '-ct', '-ci', '-cp'):
            with self.subTest(option=option):
                value = test.command_args[test.command_args.index(option) + 1]
                self.assertTrue((example_dir / value).is_file(), msg=f'{option} references missing file {value}')

    def test_amber_executable(self):
        test = load_manifest().get_test(25)
        self.assertEqual(test.executable, 'amber_MMPBSA')


if __name__ == '__main__':
    unittest.main()
