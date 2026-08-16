import ast
import logging
import unittest
from pathlib import Path

from GMXMMPBSA.exceptions import MMPBSA_Error
from GMXMMPBSA.utils import topology_mismatch_error


ROOT = Path(__file__).resolve().parents[1]


class Phase6MessageTest(unittest.TestCase):
    def test_topology_mismatch_includes_paths_and_counts(self):
        logging.disable(logging.CRITICAL)
        try:
            with self.assertRaises(MMPBSA_Error) as exc:
                topology_mismatch_error(
                    'receptor', '/work/REC.prmtop', '/work/REC.pdb', ('atoms', 1234, 1230)
                )
        finally:
            logging.disable(logging.NOTSET)

        message = str(exc.exception)
        self.assertIn('Receptor topology "/work/REC.prmtop" contains 1234 atoms', message)
        self.assertIn('receptor structure "/work/REC.pdb" contains 1230 atoms', message)
        self.assertIn('topology, structure, and index', message)

    def test_priority_messages_use_current_terminology(self):
        source = '\n'.join((ROOT / path).read_text() for path in (
            'GMXMMPBSA/main.py', 'GMXMMPBSA/make_top.py', 'GMXMMPBSA/make_top_amber.py'
        ))
        tree = ast.parse(source)
        calls = [
            ast.get_source_segment(source, node) or ''
            for node in ast.walk(tree)
            if isinstance(node, ast.Call)
        ]
        rendered = '\n'.join(calls)

        self.assertIn('prepared by cpptraj', rendered)
        self.assertIn('MPI ranks', rendered)
        self.assertIn('The mutant must', rendered)
        self.assertIn('is being used', rendered)
        self.assertNotIn('useper snapshot', rendered)
        self.assertNotIn('is been used', rendered)


if __name__ == '__main__':
    unittest.main()
