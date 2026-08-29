import csv
import tempfile
import unittest
from pathlib import Path

from GMXMMPBSA.exceptions import InputError, MMPBSA_Error
from GMXMMPBSA.membrane import (AUTOMATIC, calculate_parameters, diagnostic_paths,
                                needs_automatic_parameters, parse_atom_names,
                                read_extracted_coordinates, write_diagnostics)


class MembraneParametersTest(unittest.TestCase):
    def test_pooled_center_and_thickness_match_leaflet_means(self):
        center, thickness, diagnostics = calculate_parameters(
            {1: [-20.0, -18.0, 18.0, 20.0], 2: [-21.0, -19.0, 19.0, 21.0]},
            AUTOMATIC,
            AUTOMATIC,
        )

        self.assertAlmostEqual(center, 0.0)
        self.assertAlmostEqual(thickness, 39.0)
        self.assertEqual([row['frame'] for row in diagnostics], [1, 2])

    def test_center_and_thickness_can_be_mixed_with_manual_values(self):
        center, thickness, _ = calculate_parameters(
            {1: [-20.0, -18.0, 18.0, 20.0]},
            2.0,
            AUTOMATIC,
        )

        self.assertEqual(center, 2.0)
        self.assertAlmostEqual(thickness, 38.0)

        center, thickness, _ = calculate_parameters(
            {1: [-20.0, -18.0, 18.0, 20.0]},
            AUTOMATIC,
            42.0,
        )
        self.assertAlmostEqual(center, 0.0)
        self.assertEqual(thickness, 42.0)

    def test_atom_names_are_simple_and_semicolon_separated(self):
        self.assertEqual(parse_atom_names('P; N; P'), ('P', 'N'))
        with self.assertRaises(InputError):
            parse_atom_names(':POPC@P')

    def test_automatic_detection_requires_membrane_and_automatic_setting(self):
        self.assertTrue(needs_automatic_parameters({
            'pb': {'memopt': 1, 'mctrdz': AUTOMATIC, 'mthick': 40.0}
        }))
        self.assertFalse(needs_automatic_parameters({
            'pb': {'memopt': 0, 'mctrdz': AUTOMATIC, 'mthick': AUTOMATIC}
        }))

    def test_extracted_coordinates_pool_multiple_atom_names(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            root = Path(tmpdir)
            for atom, values in {'P': (20.0, 21.0), 'N': (-20.0, -21.0)}.items():
                for frame, z in enumerate(values, start=1):
                    (root / f'{atom}.pdb.{frame}').write_text(
                        f'HETATM    1  {atom:<3} POPC A   1    {1.0:8.3f}{2.0:8.3f}{z:8.3f}  1.00  0.00           P\n'
                    )

            coordinates = read_extracted_coordinates(root, ('P', 'N'))

        self.assertEqual(coordinates, {1: [20.0, -20.0], 2: [21.0, -21.0]})

    def test_diagnostics_are_written_to_csv_and_png(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            csv_path = Path(tmpdir) / 'membrane.csv'
            png_path = Path(tmpdir) / 'membrane.png'
            frame_coordinates = {1: [-20.0, -19.0, 19.0, 20.0], 2: [-21.0, -20.0, 20.0, 21.0]}
            _, _, diagnostics = calculate_parameters(
                frame_coordinates, AUTOMATIC, AUTOMATIC
            )
            write_diagnostics(
                diagnostics, csv_path, png_path, ('P',), AUTOMATIC, AUTOMATIC,
                0.0, 40.0, frame_coordinates=frame_coordinates,
            )

            self.assertTrue(csv_path.exists())
            self.assertTrue(png_path.exists())
            self.assertGreater(png_path.stat().st_size, 10_000)
            with csv_path.open(newline='') as handle:
                rows = list(csv.DictReader(line for line in handle if not line.startswith('#')))
            self.assertEqual(len(rows), 2)
            self.assertEqual(rows[0]['frame'], '1')

    def test_thickness_requires_both_leaflets(self):
        with self.assertRaises(MMPBSA_Error):
            calculate_parameters({1: [1.0, 1.0]}, AUTOMATIC, AUTOMATIC)

    def test_diagnostic_paths_do_not_use_the_temporary_file_prefix(self):
        csv_path, png_path = diagnostic_paths('_GMXMMPBSA_')
        self.assertEqual(csv_path.name, 'GMXMMPBSA_membrane_parameters.csv')
        self.assertEqual(png_path.name, 'GMXMMPBSA_membrane_parameters.png')


if __name__ == '__main__':
    unittest.main()
