import unittest
from inchi.determine_levels_id import InChi

class TestInChI(unittest.TestCase):
    def test_complete_identity(self):
        inchi = "InChI=1S/C2H6O/c1-2-3/h3H,2H2,1H3"
        self.assertTrue(InChi.isCompleteIdentity(inchi, inchi))

    def test_are_equal_diluted_salts(self):
        inchi_a = "InChI=1S/C2H6O/c1-2-3/h3H,2H2,1H3"
        inchi_b = "InChI=1S/C2H6O/c1-2-3/h3H,2H2,1H3"
        self.assertTrue(InChi.areEqualDisolvedSalts(inchi_a, inchi_b))

    def test_equal_salts_anion(self):
        inchi1 = "InChI=1S/C23H45NO4.ClH/c1-5-6-7-8-9-10-11-12-13-14-15-16-17-18-23(27)28-21(19-22(25)26)20-24(2,3)4;/h21H,5-20H2,1-4H3;1H"
        inchi2 = "InChI=1S/C23H45NO4.BrH/c1-5-6-7-8-9-10-11-12-13-14-15-16-17-18-23(27)28-21(19-22(25)26)20-24(2,3)4;/h21H,5-20H2,1-4H3;1H"
        self.assertTrue(InChi.areEqualDisolvedSalts(inchi1, inchi2))

    def test_equal_salts_cation(self):
        inchi1= "InChI=1S/C2H4O2.Na/c1-2(3)4;/h1H3,(H,3,4);/q;+1/p-1"
        inchi2="InChI=1S/C2H4O2.K/c1-2(3)4;/h1H3,(H,3,4);/q;+1/p-1"
        self.assertTrue(InChi.areEqualDisolvedSalts(inchi1, inchi2))

    def test_equals_without_charges(self):
        inchi1="InChI=1S/C10H16N5O13P3/c11-8-5-9(13-2-12-8)15(3-14-5)10-7(17)6(16)4(26-10)1-25-30(21,22)28-31(23,24)27-29(18,19)20/h2-4,6-7,10,16-17H,1H2,(H,21,22)(H,23,24)(H2,11,12,13)(H2,18,19,20)/t4-,6-,7-,10-/m1/s1"
        inchi2="InChI=1S/C10H16N5O13P3/c11-8-5-9(13-2-12-8)15(3-14-5)10-7(17)6(16)4(26-10)1-25-30(21,22)28-31(23,24)27-29(18,19)20/h2-4,6-7,10,16-17H,1H2,(H,21,22)(H,23,24)(H2,11,12,13)(H2,18,19,20)/p-1/t4-,6-,7-,10-/m1/s1"
        inchi3="InChI=1S/C10H16N5O13P3/c11-8-5-9(13-2-12-8)15(3-14-5)10-7(17)6(16)4(26-10)1-25-30(21,22)28-31(23,24)27-29(18,19)20/h2-4,6-7,10,16-17H,1H2,(H,21,22)(H,23,24)(H2,11,12,13)(H2,18,19,20)/p-2/t4-,6-,7-,10-/m1/s1"
        inchi4="InChI=1S/C10H16N5O13P3/c11-8-5-9(13-2-12-8)15(3-14-5)10-7(17)6(16)4(26-10)1-25-30(21,22)28-31(23,24)27-29(18,19)20/h2-4,6-7,10,16-17H,1H2,(H,21,22)(H,23,24)(H2,11,12,13)(H2,18,19,20)/p-1/t4-,6-,7-,10-/m1/s1"
        inchi5="InChI=1S/C10H16N5O13P3/c11-8-5-9(13-2-12-8)15(3-14-5)10-7(17)6(16)4(26-10)1-25-30(21,22)28-31(23,24)27-29(18,19)20/h2-4,6-7,10,16-17H,1H2,(H,21,22)(H,23,24)(H2,11,12,13)(H2,18,19,20)/p+1/t4-,6-,7-,10-/m1/s1"
        self.assertTrue(InChi.areEqualNoCharges(inchi1, inchi2))
        self.assertTrue(InChi.areEqualNoCharges(inchi3, inchi4))
        self.assertTrue(InChi.areEqualNoCharges(inchi1, inchi5))

    def test_equals_without_charges_permanent_and_optional(self):
        inchi1="InChI=1S/C5H11NO2/c1-6(2,3)4-5(7)8/h4H2,1-3H3/p+1"
        inchi2="InChI=1S/C5H11NO2/c1-6(2,3)4-5(7)8/h4H2,1-3H3"
        self.assertTrue(InChi.areEqualNoCharges(inchi1, inchi2))

    def test_are_equal_no_position_double_bond(self):
        inchi_a = "InChI=1S/C2H2/c1-2/h1-2H"  
        inchi_b = "InChI=1S/C2H2/c1-2/h1-2H"  
        self.assertTrue(InChi.areEqualNoPositionDoubleBond(inchi_a, inchi_b))

    def test_are_not_equal_no_position_double_bond(self):
        inchi_a = "InChI=1S/C2H2/c1-2/h1-2H"   
        inchi_b = "InChI=1S/CH4/h1H4"          
        self.assertFalse(InChi.areEqualNoPositionDoubleBond(inchi_a, inchi_b))

    def test_are_equal_no_isotopes(self):
        inchi_normal = "InChI=1S/C2H6O/c1-2-3/h3H,2H2,1H3"
        inchi_isotope = "InChI=1S/C2H6O/c1-2-3/h3H,2H2,1H3/i1+1"  
        self.assertTrue(InChi.areEqualNoIsotopes(inchi_normal, inchi_isotope))

    def test_are_not_equal_no_isotopes(self):
        inchi_a = "InChI=1S/C2H6O/c1-2-3/h3H,2H2,1H3"
        inchi_b = "InChI=1S/CH4O/c1-2/h2H,1H3"
        self.assertFalse(InChi.areEqualNoIsotopes(inchi_a, inchi_b))

if __name__ == "__main__":
    unittest.main()
