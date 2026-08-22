import tempfile
import unittest
from pathlib import Path
from unittest.mock import ANY, patch

from GMXMMPBSA.calculation import Calculation, EnergyCalculation, parse_qmmm_diagnostics
from GMXMMPBSA.exceptions import CalcError


class QMMMConvergenceTest(unittest.TestCase):
    def _calculation(self, input_file, output_file):
        calculation = EnergyCalculation(
            'sander', 'COM.prmtop', 'COM.inpcrd', None,
            str(input_file), str(output_file), None,
        )
        calculation.command_args = ['sander', '-o', str(output_file)]
        return calculation

    def test_unconverged_qmmm_scf_is_rejected(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            tmpdir = Path(tmpdir)
            input_file = tmpdir / 'qmmm.mdin'
            output_file = tmpdir / 'qmmm.mdout'
            input_file.write_text('&qmmm\n/\n')
            output_file.write_text(
                'QMMM: No convergence in SCF after      1 steps.\n'
                'QMMM: Job will continue with unconverged SCF. Warning energies\n'
                'QMMM: No convergence in SCF after      2 steps.\n'
            )

            calculation = self._calculation(input_file, output_file)
            with self.assertRaisesRegex(CalcError, r'did not converge in 2 step\(s\)'):
                calculation._check_qmmm_convergence()

    def test_unconverged_dftb_scf_is_rejected(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            tmpdir = Path(tmpdir)
            input_file = tmpdir / 'qmmm.mdin'
            output_file = tmpdir / 'qmmm.mdout'
            input_file.write_text('&qmmm\n/\n')
            output_file.write_text(
                'QMMM SCC-DFTB: SCC-DFTB FOR STEP     1 DID NOT CONVERGE AFTER 70 cycles.\n'
                'QMMM SCC-DFTB: Resetting Broyden mixing.\n'
                'QMMM SCC-DFTB: SCC-DFTB FOR STEP     1 DID NOT CONVERGE AFTER 140 cycles.\n'
            )

            calculation = self._calculation(input_file, output_file)
            with self.assertRaisesRegex(CalcError, r'did not converge in 2 step\(s\)'):
                calculation._check_qmmm_convergence()

    def test_qmmm_convergence_check_runs_after_sander(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            tmpdir = Path(tmpdir)
            input_file = tmpdir / 'qmmm.mdin'
            output_file = tmpdir / 'qmmm.mdout'
            input_file.write_text('&qmmm\n/\n')
            output_file.write_text('QMMM SCF converged.\n')

            calculation = self._calculation(input_file, output_file)
            with patch.object(Calculation, 'run') as base_run:
                calculation.run(0)
            base_run.assert_called_once_with(ANY, 0, ANY, ANY)

    def test_non_qmmm_output_is_not_checked_for_qm_messages(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            tmpdir = Path(tmpdir)
            input_file = tmpdir / 'gb.mdin'
            output_file = tmpdir / 'gb.mdout'
            input_file.write_text('&cntrl\n/\n')
            output_file.write_text('No convergence in SCF after      1 steps.\n')

            calculation = self._calculation(input_file, output_file)
            calculation._check_qmmm_convergence()

    def test_missing_method_parameters_include_method_and_atomic_number(self):
        diagnostics = parse_qmmm_diagnostics(
            'QMMM: Atom number: 74 has atomic number 12.\n'
            'QMMM: There are no RM1 parameters for this element. Sorry.\n'
        )

        self.assertEqual(len(diagnostics), 1)
        self.assertEqual(diagnostics[0].code, 'qm_parameter_missing')
        self.assertIn('RM1', diagnostics[0].message)
        self.assertIn('atomic number 12', diagnostics[0].message)

    def test_missing_method_parameters_support_explicit_atomic_number(self):
        diagnostics = parse_qmmm_diagnostics(
            'QMMM: There are no PM3-MAIS parameters for atomic number 15. Sorry.\n'
        )

        self.assertEqual(len(diagnostics), 1)
        self.assertIn('PM3-MAIS', diagnostics[0].message)
        self.assertIn('atomic number 15', diagnostics[0].message)

    def test_dispersion_and_dftb_file_failures_are_classified(self):
        diagnostics = parse_qmmm_diagnostics(
            'SANDER BOMB in subroutine get_c6_req (dh_correction_module)\n'
            'Parameters for dispersion correction are not available for this atom.\n'
            'Missing file:\n'
            '/tmp/slko/mio-1-1/P-P.skf\n',
            qm_theory='PM6-DH+',
        )

        self.assertEqual(
            [diagnostic.code for diagnostic in diagnostics],
            ['dispersion_parameter_missing', 'dftb_parameter_file_missing'],
        )
        self.assertIn('PM6-DH+', diagnostics[0].message)
        self.assertIn('/tmp/slko/mio-1-1/P-P.skf', diagnostics[1].message)

    def test_numerical_derivative_message_is_a_warning(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            tmpdir = Path(tmpdir)
            input_file = tmpdir / 'qmmm.mdin'
            output_file = tmpdir / 'qmmm.mdout'
            input_file.write_text('&qmmm\nqm_theory="PM6"\n/\n')
            output_file.write_text(
                'QMMM: Analytical derivatives for d orbitals are not supported.\n'
            )

            calculation = self._calculation(input_file, output_file)
            with self.assertLogs(level='WARNING') as logs:
                calculation._check_qmmm_convergence()
            self.assertIn('numerical derivatives', logs.output[0])

    def test_sander_failure_is_replaced_with_qmmm_diagnostic(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            tmpdir = Path(tmpdir)
            input_file = tmpdir / 'qmmm.mdin'
            output_file = tmpdir / 'qmmm.mdout'
            input_file.write_text('&qmmm\nqm_theory="RM1"\n/\n')
            output_file.write_text(
                'QMMM: Atom number: 74 has atomic number 12.\n'
                'QMMM: There are no RM1 parameters for this element. Sorry.\n'
            )

            calculation = self._calculation(input_file, output_file)
            with patch.object(Calculation, 'run', side_effect=CalcError('sander failed')):
                with self.assertRaisesRegex(CalcError, r'qm_parameter_missing.*RM1.*atomic number 12'):
                    calculation.run(0)


if __name__ == '__main__':
    unittest.main()
