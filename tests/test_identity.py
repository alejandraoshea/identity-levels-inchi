import unittest
from inchi.determine_levels_id import isCompleteIdentity, areEqualDisolvedSalts, areEqualNoPositionDoubleBond, areEqualNoIsotopes

class TestInChI(unittest.TestCase):
    def test_complete_identity(self):
        inchi = "InChI=1S/C2H6O/c1-2-3/h3H,2H2,1H3"
        self.assertTrue(isCompleteIdentity(inchi, inchi))

    def test_are_equal_diluted_salts(self):
        inchi_a = "InChI=1S/C2H6O/c1-2-3/h3H,2H2,1H3"
        inchi_b = "InChI=1S/C2H6O/c1-2-3/h3H,2H2,1H3"
        self.assertTrue(areEqualDisolvedSalts(inchi_a, inchi_b))

    def test_are_equal_no_position_double_bond(self):
        inchi_a = "InChI=1S/C2H2/c1-2/h1-2H"  
        inchi_b = "InChI=1S/C2H2/c1-2/h1-2H"  
        self.assertTrue(areEqualNoPositionDoubleBond(inchi_a, inchi_b))

    def test_are_not_equal_no_position_double_bond(self):
        inchi_a = "InChI=1S/C2H2/c1-2/h1-2H"   
        inchi_b = "InChI=1S/CH4/h1H4"          
        self.assertFalse(areEqualNoPositionDoubleBond(inchi_a, inchi_b))

    def test_are_equal_no_isotopes(self):
        inchi_normal = "InChI=1S/C2H6O/c1-2-3/h3H,2H2,1H3"
        inchi_isotope = "InChI=1S/C2H6O/c1-2-3/h3H,2H2,1H3/i1+1"  
        self.assertTrue(areEqualNoIsotopes(inchi_normal, inchi_isotope))

    def test_are_not_equal_no_isotopes(self):
        inchi_a = "InChI=1S/C2H6O/c1-2-3/h3H,2H2,1H3"
        inchi_b = "InChI=1S/CH4O/c1-2/h2H,1H3"
        self.assertFalse(areEqualNoIsotopes(inchi_a, inchi_b))

if __name__ == "__main__":
    unittest.main()
