import unittest

from GMXMMPBSA.make_top import CheckMakeTop


class GromacsRadiusCompatibilityTest(unittest.TestCase):
    @staticmethod
    def _checker(pbradii, igb=8, gbrun=True):
        checker = object.__new__(CheckMakeTop)
        checker.INPUT = {
            'general': {'PBRadii': pbradii},
            'gb': {'gbrun': gbrun, 'igb': igb},
        }
        return checker

    def test_warns_for_nonconventional_gb_radius_pairing(self):
        checker = self._checker(pbradii=3, igb=8)

        with self.assertLogs(level='WARNING') as messages:
            checker._warn_gmx_gb_radius_compatibility()

        self.assertIn("PBRadii='mbondi2' is selected", messages.output[0])
        self.assertIn("igb=8", messages.output[0])
        self.assertIn("'mbondi3' radii", messages.output[0])

    def test_does_not_warn_for_conventional_gb_radius_pairing(self):
        checker = self._checker(pbradii=4, igb=8)

        with self.assertNoLogs(level='WARNING'):
            checker._warn_gmx_gb_radius_compatibility()

    def test_does_not_warn_when_gb_is_not_selected(self):
        checker = self._checker(pbradii=3, igb=8, gbrun=False)

        with self.assertNoLogs(level='WARNING'):
            checker._warn_gmx_gb_radius_compatibility()


if __name__ == '__main__':
    unittest.main()
