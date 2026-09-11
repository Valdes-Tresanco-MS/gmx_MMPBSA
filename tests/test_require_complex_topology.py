import unittest
from types import SimpleNamespace
from unittest.mock import patch

from GMXMMPBSA import make_top
from GMXMMPBSA.exceptions import MMPBSA_Error


class RequireComplexTopologyTest(unittest.TestCase):
    def _checker(self, **file_overrides):
        files = SimpleNamespace(
            complex_tpr='com.tpr',
            complex_index='index.ndx',
            complex_trajs=['com.xtc'],
            complex_groups=['Protein', 'LIG'],
            complex_top=None,
            receptor_tpr=None,
            receptor_trajs=None,
            receptor_top=None,
            ligand_tpr=None,
            ligand_trajs=None,
            ligand_top=None,
            prefix='_GMXMMPBSA_',
        )
        for key, value in file_overrides.items():
            setattr(files, key, value)
        checker = object.__new__(make_top.CheckMakeTop)
        checker.FILES = files
        return checker

    def test_missing_complex_top_raises(self):
        checker = self._checker()
        with self.assertRaises(MMPBSA_Error) as raised:
            checker.checkFiles()
        self.assertIn('-cp', str(raised.exception))
        self.assertIn('tleap', str(raised.exception).lower())

    def test_complex_top_alone_is_accepted(self):
        checker = self._checker(complex_top='topol.top')
        checker.checkFiles()

    def test_mt_receptor_requires_receptor_top(self):
        checker = self._checker(complex_top='topol.top', receptor_trajs=['rec.xtc'])
        with self.assertRaises(MMPBSA_Error) as raised:
            checker.checkFiles()
        self.assertIn('-rp', str(raised.exception))

    def test_mt_ligand_requires_ligand_top(self):
        checker = self._checker(complex_top='topol.top', ligand_tpr='lig.tpr')
        with self.assertRaises(MMPBSA_Error) as raised:
            checker.checkFiles()
        self.assertIn('-lp', str(raised.exception))

    def test_build_topology_always_uses_gmxtop2prmtop(self):
        checker = self._checker(complex_top='topol.top')
        checker.explicit_waters_mask = ''
        checker.INPUT = {
            'decomp': {'decomprun': False},
            'gb': {'ifqnt': 0, 'com_qmmask': ''},
        }
        checker.gmx2pdb = lambda: None
        checker._resolve_explicit_waters_mask = lambda: None
        checker.cleanup_trajs = lambda: None
        checker.gmxtop2prmtop = lambda: ('COM.prmtop', 'REC.prmtop', 'LIG.prmtop', None, None, None)

        with patch.object(checker, 'pdb2prmtop') as pdb2prmtop, \
                patch.object(checker, 'makeToptleap') as make_tleap:
            tops = checker.buildTopology()

        self.assertEqual(tops[0], 'COM.prmtop')
        pdb2prmtop.assert_not_called()
        make_tleap.assert_not_called()


if __name__ == '__main__':
    unittest.main()
