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
    validate_output_paths,
)
from GMXMMPBSA.exceptions import MMPBSA_Error
from GMXMMPBSA.utils import _get_dup_args, get_index_groups


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

    def test_accepts_legacy_membrane_test_selector(self):
        parser = testparser.parse_args(['-t', '11'])
        self.assertEqual(parser.test, ['11'])

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


class CsvOutputFilenameTest(unittest.TestCase):
    def test_csv_outputs_follow_summary_names_by_default(self):
        args = parser.parse_args(['-o', 'results.dat', '-do', 'decomposition.dat'])

        self.assertEqual(args.energyout, 'results.csv')
        self.assertEqual(args.dec_energies, 'decomposition.csv')

    def test_explicit_csv_names_override_the_defaults(self):
        args = parser.parse_args([
            '-o', 'results.dat', '-do', 'decomposition.dat',
            '-eo', 'energies_for_analysis.csv', '-deo', 'residue_energies.csv',
        ])

        self.assertEqual(args.energyout, 'energies_for_analysis.csv')
        self.assertEqual(args.dec_energies, 'residue_energies.csv')

    def test_defaulting_helper_supports_amber_parser_options(self):
        args = amber_parser.parse_args(['-o', 'amber.out', '-do', 'amber_decomp.out'])

        self.assertEqual(args.energyout, 'amber.csv')
        self.assertEqual(args.dec_energies, 'amber_decomp.csv')

    def test_csv_summary_uses_distinct_automatic_name_for_both_engines(self):
        for cli in (parser, amber_parser):
            args = cli.parse_args(['-o', 'results.csv', '-do', 'residues.csv'])
            self.assertEqual(args.energyout, 'results.frames.csv')
            self.assertEqual(args.dec_energies, 'residues.frames.csv')
            validate_output_paths(args, True)

    def test_active_output_collisions_leave_existing_file_untouched(self):
        for cli in (parser, amber_parser):
            with TemporaryDirectory() as directory:
                dest = Path(directory) / 'summary.dat'
                dest.write_text('preserve me')
                alias = Path(directory) / 'alias.dat'
                alias.symlink_to(dest)
                for option, value in [('-eo', str(dest.parent / '.' / dest.name)),
                                      ('-eo', str(alias)), ('-do', str(dest)), ('-deo', str(dest))]:
                    args = cli.parse_args(['-o', str(dest), option, value])
                    with self.assertRaisesRegex(MMPBSA_Error, 'Output paths for'):
                        validate_output_paths(args, True)
                    self.assertEqual(dest.read_text(), 'preserve me')
                args = cli.parse_args(['-o', str(dest), '-do', str(dest)])
                validate_output_paths(args, False)

    def test_writer_rejects_collision_before_truncating(self):
        from types import SimpleNamespace
        from GMXMMPBSA.output_file import write_outputs
        with TemporaryDirectory() as directory:
            path = Path(directory) / 'results.csv'
            path.write_text('preserve me')
            files = parser.parse_args(['-o', str(path), '-eo', str(path)])
            app = SimpleNamespace(FILES=files, INPUT={'decomp': {'decomprun': False}},
                                  normal_system=None, mut_str='', stability=False)
            with self.assertRaisesRegex(MMPBSA_Error, 'Output paths for -o and -eo'):
                write_outputs(app)
            self.assertEqual(path.read_text(), 'preserve me')


class AmberComplexStructureOptionTest(unittest.TestCase):
    def test_complex_structure_is_not_required_or_advertised(self):
        self.assertFalse(hasattr(amber_parser.parse_args([]), 'complex_str'))
        self.assertNotIn('-cs', amber_parser.format_help())

    def test_complex_structure_option_is_not_part_of_amber_interface(self):
        with self.assertRaises(MMPBSA_Error) as exc:
            amber_parser.parse_args(['-cs', 'complex.inpcrd'])

        self.assertIn('unrecognized arguments: -cs', str(exc.exception))


class GromacsGroupArgumentTest(unittest.TestCase):
    def test_receptor_and_ligand_groups_accept_numbers(self):
        args = parser.parse_args(['-rg', '1', '-lg', '2'])

        self.assertEqual(args.receptor_group, 1)
        self.assertEqual(args.ligand_group, 2)

    def test_receptor_and_ligand_groups_accept_names(self):
        args = parser.parse_args(['-rg', 'Protein chain A', '-lg', 'Ligand'])

        self.assertEqual(args.receptor_group, 'Protein chain A')
        self.assertEqual(args.ligand_group, 'Ligand')

    def test_named_receptor_and_ligand_groups_are_resolved(self):
        with TemporaryDirectory() as tmpdir:
            index_file = Path(tmpdir) / 'index.ndx'
            index_file.write_text('[ System ]\n1 2 3\n[ Protein chain A ]\n1 2\n[ Ligand ]\n3\n')

            self.assertEqual(get_index_groups(index_file, 'Protein chain A'), (1, 'Protein chain A'))
            self.assertEqual(get_index_groups(index_file, 'Ligand'), (2, 'Ligand'))

    def test_numeric_receptor_and_ligand_groups_are_resolved(self):
        with TemporaryDirectory() as tmpdir:
            index_file = Path(tmpdir) / 'index.ndx'
            index_file.write_text('[ System ]\n1 2 3\n[ Protein ]\n1 2\n[ Ligand ]\n3\n')

            self.assertEqual(get_index_groups(index_file, 1), (1, 'Protein'))
            self.assertEqual(get_index_groups(index_file, 2), (2, 'Ligand'))

    def test_duplicate_group_names_require_a_numeric_selection(self):
        with TemporaryDirectory() as tmpdir:
            index_file = Path(tmpdir) / 'index.ndx'
            index_file.write_text('[ System ]\n1 2 3\n[ Ligand ]\n2\n[ Ligand ]\n3\n')

            logging.disable(logging.CRITICAL)
            try:
                with self.assertRaisesRegex(
                    MMPBSA_Error,
                    r"Index group name 'Ligand' is ambiguous because it occurs at group numbers 1, 2\. "
                    r'Select the group by number instead\.',
                ):
                    get_index_groups(index_file, 'Ligand')
            finally:
                logging.disable(logging.NOTSET)

            self.assertEqual(get_index_groups(index_file, 2), (2, 'Ligand'))

    def test_unknown_group_name_and_out_of_range_number_are_rejected(self):
        with TemporaryDirectory() as tmpdir:
            index_file = Path(tmpdir) / 'index.ndx'
            index_file.write_text('[ Protein ]\n1 2\n')

            logging.disable(logging.CRITICAL)
            try:
                for group in ('Missing', -1, 1):
                    with self.subTest(group=group):
                        with self.assertRaisesRegex(MMPBSA_Error, 'Define a valid index group'):
                            get_index_groups(index_file, group)
            finally:
                logging.disable(logging.NOTSET)


class AmberMultipleTrajectoryArgumentTest(unittest.TestCase):
    def test_reused_amber_masks_are_not_reported_as_duplicate_arguments(self):
        _get_dup_args([
            '-cm', ':1-166', ':167-242',
            '-rm', ':1-166',
            '-lm', ':167-242',
        ])
