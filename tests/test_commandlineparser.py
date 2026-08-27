import logging
import unittest
from pathlib import Path
from tempfile import TemporaryDirectory

from GMXMMPBSA.commandlineparser import (
    AMBER_TRAJECTORY_FORMATS,
    amber_parser,
    amber_trajectory,
    parser,
    testparser,
)
from GMXMMPBSA.exceptions import MMPBSA_Error
from GMXMMPBSA.utils import _get_dup_args


class AmberTrajectoryTypeTest(unittest.TestCase):
    def test_accepts_cpptraj_readable_trajectory_formats(self):
        with TemporaryDirectory() as tmpdir:
            root = Path(tmpdir)
            for suffix in AMBER_TRAJECTORY_FORMATS:
                traj = root / f'traj{suffix}'
                traj.touch()

                self.assertEqual(amber_trajectory(traj.as_posix()), traj.as_posix())

    def test_rejects_unsupported_trajectory_format(self):
        with TemporaryDirectory() as tmpdir:
            traj = Path(tmpdir) / 'traj.xyz'
            traj.touch()

            logging.disable(logging.CRITICAL)
            try:
                with self.assertRaises(MMPBSA_Error):
                    amber_trajectory(traj.as_posix())
            finally:
                logging.disable(logging.NOTSET)


class TestParserSelectorTest(unittest.TestCase):
    def test_accepts_named_test_selector(self):
        parser = testparser.parse_args(['-t', 'gbnsr6'])
        self.assertEqual(parser.test, ['gbnsr6'])

    def test_rejects_unknown_test_selector(self):
        with self.assertRaises(SystemExit) as exc:
            testparser.parse_args(['-t', 'not-a-test'])
        self.assertEqual(exc.exception.code, 2)


class ProgressStyleParserTest(unittest.TestCase):
    def test_progress_style_defaults_to_auto(self):
        self.assertEqual(parser.parse_args([]).progress_style, 'auto')
        self.assertEqual(amber_parser.parse_args([]).progress_style, 'auto')

    def test_progress_style_can_be_selected(self):
        self.assertEqual(parser.parse_args(['--progress-style', 'classic']).progress_style, 'classic')
        self.assertEqual(amber_parser.parse_args(['--progress-style', 'none']).progress_style, 'none')


class AmberComplexStructureOptionTest(unittest.TestCase):
    def test_complex_structure_is_not_required_or_advertised(self):
        self.assertFalse(hasattr(amber_parser.parse_args([]), 'complex_str'))
        self.assertNotIn('-cs', amber_parser.format_help())

    def test_complex_structure_option_is_not_part_of_amber_interface(self):
        with self.assertRaises(MMPBSA_Error) as exc:
            amber_parser.parse_args(['-cs', 'complex.inpcrd'])

        self.assertIn('unrecognized arguments: -cs', str(exc.exception))


class AmberMultipleTrajectoryArgumentTest(unittest.TestCase):
    def test_reused_amber_masks_are_not_reported_as_duplicate_arguments(self):
        _get_dup_args([
            '-cm', ':1-166', ':167-242',
            '-rm', ':1-166',
            '-lm', ':167-242',
        ])
