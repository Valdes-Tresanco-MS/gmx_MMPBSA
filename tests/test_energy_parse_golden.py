"""Fixture-based GB energy parse golden (kcal/mol component means).

This is not a full AmberTools end-to-end regression. It locks the
mdout → EnergyVector → BindingStatistics path so GB totals cannot drift
silently. See scripts/validation/README.md for adding Amber binaries goldens.
"""

import tempfile
import unittest
from pathlib import Path

from GMXMMPBSA.amber_outputs import BindingStatistics, GBout


def _gb_input():
    return {
        'general': {'temperature': 298.15, 'startframe': 1, 'interval': 1},
        'gb': {'surften': 1.0, 'surfoff': 0.0},
        'pb': {'sander_apbs': 0},
        'nmode': {'nmstartframe': 1, 'nminterval': 1},
    }


def _write_gb_mdout(path, frames):
    """Write a minimal sander-style GB energy block sequence."""
    lines = []
    for bond, angle, dihed, vdw, eel, egb, vdw14, eel14 in frames:
        lines.extend(
            [
                f' BOND    = {bond:9.4f}  ANGLE   = {angle:9.4f}  DIHED      = {dihed:9.4f}',
                f' VDWAALS = {vdw:9.4f}  EEL     = {eel:9.4f}  EGB        = {egb:9.4f}',
                f' 1-4 VDW = {vdw14:9.4f}  1-4 EEL = {eel14:9.4f}',
                '',
            ]
        )
    Path(path).write_text('\n'.join(lines) + '\n')


def _write_surf(path, areas):
    Path(path).write_text('#Frame Area\n' + ''.join(f'{i} {a}\n' for i, a in enumerate(areas, 1)))


class GbParseGoldenTest(unittest.TestCase):
    def test_delta_component_means_match_hand_calculated_kcal(self):
        # Two STP frames; bonded terms cancel in the delta for each frame.
        # Frame averages (kcal/mol):
        #   ΔVDWAALS=-30, ΔEEL=-50, ΔEGB=-15, ΔESURF=0.5 → ΔGGAS=-80, ΔGSOLV=-14.5, ΔTOTAL=-94.5
        com_frames = [
            (15.0, 25.0, 35.0, -100.0, -200.0, -50.0, 5.0, 6.0),
            (15.0, 25.0, 35.0, -102.0, -198.0, -52.0, 5.0, 6.0),
        ]
        rec_frames = [
            (10.0, 20.0, 30.0, -40.0, -80.0, -20.0, 3.0, 4.0),
            (10.0, 20.0, 30.0, -41.0, -79.0, -21.0, 3.0, 4.0),
        ]
        lig_frames = [
            (5.0, 5.0, 5.0, -30.0, -70.0, -15.0, 2.0, 2.0),
            (5.0, 5.0, 5.0, -31.0, -69.0, -16.0, 2.0, 2.0),
        ]
        # surften=1.0 → ESURF equals the printed "area" column.
        com_surf = [2.0, 2.0]
        rec_surf = [1.0, 1.0]
        lig_surf = [0.5, 0.5]

        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            for mol, frames, surf in (
                ('complex', com_frames, com_surf),
                ('receptor', rec_frames, rec_surf),
                ('ligand', lig_frames, lig_surf),
            ):
                _write_gb_mdout(root / f'{mol}_gb.mdout.0', frames)
                _write_surf(root / f'{mol}_gb_surf.dat.0', surf)

            INPUT = _gb_input()
            com = GBout('complex', INPUT)
            rec = GBout('receptor', INPUT)
            lig = GBout('ligand', INPUT)
            com.parse_from_file(str(root / 'complex_gb.mdout'), num_files=1, numframes=2)
            rec.parse_from_file(str(root / 'receptor_gb.mdout'), num_files=1, numframes=2)
            lig.parse_from_file(str(root / 'ligand_gb.mdout'), num_files=1, numframes=2)

            delta = BindingStatistics(com, rec, lig, traj_protocol='STP')

        self.assertAlmostEqual(float(delta['VDWAALS'].mean()), -30.0, places=6)
        self.assertAlmostEqual(float(delta['EEL'].mean()), -50.0, places=6)
        self.assertAlmostEqual(float(delta['EGB'].mean()), -15.0, places=6)
        self.assertAlmostEqual(float(delta['ESURF'].mean()), 0.5, places=6)
        self.assertAlmostEqual(float(delta['BOND'].mean()), 0.0, places=6)
        self.assertAlmostEqual(float(delta['GGAS'].mean()), -80.0, places=6)
        self.assertAlmostEqual(float(delta['GSOLV'].mean()), -14.5, places=6)
        self.assertAlmostEqual(float(delta['TOTAL'].mean()), -94.5, places=6)


if __name__ == '__main__':
    unittest.main()
