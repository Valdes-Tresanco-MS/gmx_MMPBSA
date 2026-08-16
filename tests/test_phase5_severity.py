import ast
import unittest
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]


def _logging_levels(path, fragment):
    source = path.read_text()
    tree = ast.parse(source)
    levels = []
    for node in ast.walk(tree):
        if not isinstance(node, ast.Call) or not isinstance(node.func, ast.Attribute):
            continue
        if not (isinstance(node.func.value, ast.Name) and node.func.value.id == 'logging'):
            continue
        call_text = ast.get_source_segment(source, node) or ''
        if fragment in call_text:
            levels.append(node.func.attr)
    return levels


class ExpectedActionSeverityTest(unittest.TestCase):
    def test_automatic_actions_are_info(self):
        expected_info = {
            'GMXMMPBSA/createinput.py': [
                'Setting complex arad', 'Setting receptor arad', 'Setting ligand arad',
            ],
            'GMXMMPBSA/make_top.py': [
                'Setting qmcharge_', 'Using user-defined qmcharge_',
                'Generating a receptor file internally',
                'Stability calculation mode does not need',
                'Assigning missing chain IDs',
            ],
            'GMXMMPBSA/make_top_amber.py': [
                'Setting qmcharge_', 'Using user-defined qmcharge_',
                'Generating a receptor file internally',
                'Stability calculation mode does not need',
                'Assigning missing chain IDs',
            ],
            'GMXMMPBSA/main.py': ['Preparing GBNSR6 topology copies'],
        }

        for relative_path, fragments in expected_info.items():
            path = ROOT / relative_path
            for fragment in fragments:
                with self.subTest(path=relative_path, fragment=fragment):
                    self.assertTrue(_logging_levels(path, fragment))
                    self.assertEqual(_logging_levels(path, fragment), ['info'] * len(_logging_levels(path, fragment)))

    def test_scientific_or_fallback_consequences_remain_warnings(self):
        expected_warning = {
            'GMXMMPBSA/make_top.py': [
                'invalid DIHEDRAL_PERIODICITY',
                'Could not identify GB radii',
                'EXPLICIT_WATERS_EXTRA_POINTS',
                'omits CMAP energy terms',
                'Reassigning existing chain IDs',
            ],
            'GMXMMPBSA/make_top_amber.py': [
                'invalid DIHEDRAL_PERIODICITY',
                'Could not identify GB radii',
                'Reassigning existing chain IDs',
            ],
        }

        for relative_path, fragments in expected_warning.items():
            path = ROOT / relative_path
            for fragment in fragments:
                with self.subTest(path=relative_path, fragment=fragment):
                    self.assertEqual(_logging_levels(path, fragment), ['warning'] * len(_logging_levels(path, fragment)))


if __name__ == '__main__':
    unittest.main()
