import unittest
from tempfile import TemporaryDirectory
from pathlib import Path

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
