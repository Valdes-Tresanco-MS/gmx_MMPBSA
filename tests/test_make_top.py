import unittest
from tempfile import TemporaryDirectory
from pathlib import Path
from types import SimpleNamespace
from unittest.mock import patch

from GMXMMPBSA import make_top
from GMXMMPBSA.topology_preprocess import GromacsTopologyPreprocessor, comment_gromacs_cmap


class CommentGromacsCmapTest(unittest.TestCase):
    def _filter(self, lines):
        in_cmap = False
        output = []
        for line in lines:
            line, in_cmap, _ = comment_gromacs_cmap(line, in_cmap)
            output.append(line)
        return output

    def test_comments_complete_cmap_block(self):
        lines = [
            '[ dihedrals ]\n',
            '1 2 3 4 1\n',
            '[ cmap ]\n',
            '; ai aj ak al am funct\n',
            '1 2 3 4 5 1\n',
            '[ molecules ]\n',
            'Protein 1\n',
        ]

        self.assertEqual(
            self._filter(lines),
            [
                '[ dihedrals ]\n',
                '1 2 3 4 1\n',
                ';[ cmap ]\n',
                '; ai aj ak al am funct\n',
                ';1 2 3 4 5 1\n',
                '[ molecules ]\n',
                'Protein 1\n',
            ],
        )

    def test_detects_whitespace_and_case_variants(self):
        lines = [
            '  [  CMAP  ]\n',
            '1 2 3 4 5 1\n',
            '[ atoms ]\n',
        ]

        self.assertEqual(
            self._filter(lines),
            [
                ';  [  CMAP  ]\n',
                ';1 2 3 4 5 1\n',
                '[ atoms ]\n',
            ],
        )

    def test_ignores_already_commented_cmap_header(self):
        lines = [
            ';[ cmap ]\n',
            '1 2 3 4 5 1\n',
            '[ atoms ]\n',
        ]

        self.assertEqual(self._filter(lines), lines)

    def test_preserves_preprocessor_directives_after_cmap(self):
        lines = [
            '[ cmap ]\n',
            '1 2 3 4 5 1\n',
            '#ifdef POSRES\n',
            '[ position_restraints ]\n',
            '1 1 1000 1000 1000\n',
            '#endif\n',
        ]

        self.assertEqual(
            self._filter(lines),
            [
                ';[ cmap ]\n',
                ';1 2 3 4 5 1\n',
                '#ifdef POSRES\n',
                '[ position_restraints ]\n',
                '1 1 1000 1000 1000\n',
                '#endif\n',
            ],
        )

    def test_preprocesses_local_includes_recursively(self):
        with TemporaryDirectory() as tmpdir:
            root = Path(tmpdir)
            top = root / 'topol.top'
            protein = root / 'topol_Protein.itp'
            posre = root / 'posre_Protein.itp'

            top.write_text(
                '#include "amber19sb.ff/forcefield.itp"\n'
                '#include "topol_Protein.itp"\n'
                '[ molecules ]\n'
                'Protein 1\n'
                'SOL 10\n'
            )
            protein.write_text(
                '[ moleculetype ]\n'
                'Protein 3\n'
                '[ cmap ]\n'
                '1 2 3 4 5 1\n'
                '[ atoms ]\n'
                '1 C 1 ALA CA 1 0 12.01\n'
                '#include "posre_Protein.itp"\n'
            )
            posre.write_text('[ position_restraints ]\n1 1 1000 1000 1000\n')

            preprocessor = GromacsTopologyPreprocessor()
            sanitized_top = preprocessor.preprocess(top, True, ['SOL'])

            sanitized_text = sanitized_top.read_text()
            self.assertIn('#include "amber19sb.ff/forcefield.itp"', sanitized_text)
            self.assertNotIn('#include "topol_Protein.itp"', sanitized_text)
            self.assertNotIn('SOL 10', sanitized_text)
            self.assertEqual(len(sanitized_text.splitlines()), 4)

            sanitized_protein = next(
                path for path in preprocessor.created_files
                if path.name.startswith('_temp_topol_Protein_')
            )
            protein_text = sanitized_protein.read_text()
            self.assertIn(';[ cmap ]', protein_text)
            self.assertIn(';1 2 3 4 5 1', protein_text)
            self.assertIn('[ atoms ]', protein_text)
            self.assertNotIn('#include "posre_Protein.itp"', protein_text)

            for temp_file in preprocessor.created_files:
                temp_file.unlink(missing_ok=True)

    def test_cmap_omission_is_reported_as_approximation_warning(self):
        preprocessor = SimpleNamespace(
            cmap_found=True,
            created_files=[],
            preprocess=lambda *args: Path('sanitized.top'),
        )
        residue = SimpleNamespace(idx=0, number=0)
        topology = SimpleNamespace(
            atoms=[SimpleNamespace(residue=residue)],
            residues=[residue],
            strip=lambda selection: None,
        )

        with patch.object(make_top, 'GromacsTopologyPreprocessor', return_value=preprocessor), \
                patch.object(make_top.parmed.gromacs, 'GromacsTopologyFile', return_value=topology):
            with self.assertLogs(level='WARNING') as logs:
                make_top.CheckMakeTop.cleantop('topol.top', [1])

        self.assertIn('omits CMAP energy terms', logs.output[0])
        self.assertIn('not an issue', logs.output[0])
        self.assertIn('MTP', logs.output[0])
        self.assertTrue(logs.output[0].startswith('WARNING:'))


class ChainAssignmentTest(unittest.TestCase):
    def test_terminal_residue_number_gap_does_not_index_past_residue_map(self):
        checker = object.__new__(make_top.CheckMakeTop)
        checker.resl = [
            make_top.Residue(1, 153, '', 'R', 1, 'GLY'),
            make_top.Residue(2, 155, '', 'L', 1, 'HEME'),
        ]
        complex_structure = SimpleNamespace(residues=[
            SimpleNamespace(chain='', number=153, atoms=[], name='GLY'),
            SimpleNamespace(chain='', number=155, atoms=[], name='HEME'),
        ])
        receptor_structure = SimpleNamespace(residues=[SimpleNamespace(chain='')])
        ligand_structure = SimpleNamespace(residues=[SimpleNamespace(chain='')])

        checker._assign_chains_IDs(complex_structure, receptor_structure, ligand_structure)

        self.assertEqual([res.chain for res in complex_structure.residues], ['A', 'B'])
        self.assertEqual(receptor_structure.residues[0].chain, 'A')
        self.assertEqual(ligand_structure.residues[0].chain, 'B')

    def test_preserves_relative_path_for_subdirectory_includes(self):
        with TemporaryDirectory() as tmpdir:
            root = Path(tmpdir)
            forcefield = root / 'amber19sb.ff'
            forcefield.mkdir()
            top = root / 'topol.top'
            ff = forcefield / 'forcefield.itp'

            top.write_text('#include "amber19sb.ff/forcefield.itp"\n')
            ff.write_text('[ cmap ]\n1 2 3 4 5 1\n')

            preprocessor = GromacsTopologyPreprocessor()
            sanitized_top = preprocessor.preprocess(top)

            self.assertIn('#include "amber19sb.ff/_temp_forcefield_', sanitized_top.read_text())

            for temp_file in preprocessor.created_files:
                temp_file.unlink(missing_ok=True)


if __name__ == '__main__':
    unittest.main()
