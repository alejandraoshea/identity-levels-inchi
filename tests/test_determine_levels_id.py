import unittest
from inchi.determine_levels_id import InChi
from inchi.inchi_layers_enum import InchiLayers


class TestIdentityLevels(unittest.TestCase):

    @staticmethod
    def check(layer: InchiLayers, inchi1: str, inchi2: str):
        """Helper: run classifier and assert only the given layer is true."""
        results = InChi.get_ids(inchi1, inchi2)

        assert layer in results, "Layer missing in result map"

        # True for the intended layer
        assert results[layer] is True, f"{layer} should be True"

        # All others must be False
        for other, value in results.items():
            if other != layer:
                assert value is False, f"{other} should be False, but got True"

    # A — Complete Identity
    def test_complete_identity(self):
        i1 = "InChI=1S/CH4/h1H4"
        i2 = "InChI=1S/CH4/h1H4"
        self.check(InchiLayers.COMPLETE_IDENTITY, i1, i2)

    # B — Identity independent of salts
    def test_independent_salts_palmitoyl_carnitine(self):
        i1 = "InChI=1S/C23H45NO4.ClH/c1-5-6-7-8-9-10-11-12-13-14-15-16-17-18-23(27)28-21(19-22(25)26)20-24(2,3)4;/h21H,5-20H2,1-4H3;1H"
        i2 = "InChI=1S/C23H45NO4.BrH/c1-5-6-7-8-9-10-11-12-13-14-15-16-17-18-23(27)28-21(19-22(25)26)20-24(2,3)4;/h21H,5-20H2,1-4H3;1H"
        self.check(InchiLayers.INDEPENDENT_SALTS, i1, i2)

    def test_independent_salts_acetate(self):
        i1 = "InChI=1S/C2H4O2.Na/c1-2(3)4;/h1H3,(H,3,4);/q;+1/p-1"
        i2 = "InChI=1S/C2H4O2.K/c1-2(3)4;/h1H3,(H,3,4);/q;+1/p-1"
        self.check(InchiLayers.INDEPENDENT_SALTS, i1, i2)

    # C — Identity independent of charge
    def test_independent_charges_atp(self):
        i1 = "InChI=1S/C10H16N5O13P3/c11-8-5-9(13-2-12-8)15(3-14-5)10-7(17)6(16)4(26-10)1-25-30(21,22)28-31(23,24)27-29(18,19)20/h2-4,6-7,10,16-17H,1H2,(H,21,22)(H,23,24)(H2,11,12,13)(H2,18,19,20)/t4-,6-,7-,10-/m1/s1"
        i2 = "InChI=1S/C10H16N5O13P3/c11-8-5-9(13-2-12-8)15(3-14-5)10-7(17)6(16)4(26-10)1-25-30(21,22)28-31(23,24)27-29(18,19)20/h2-4,6-7,10,16-17H,1H2,(H,21,22)(H,23,24)(H2,11,12,13)(H2,18,19,20)/p-1/t4-,6-,7-,10-/m1/s1"
        self.check(InchiLayers.INDEPENDENT_CHARGES, i1, i2)

    def test_independent_charges_permanent_plus_optional_minus(self):
        i1 = "InChI=1S/C5H11NO2/c1-6(2,3)4-5(7)8/h4H2,1-3H3/p+1"
        i2 = "InChI=1S/C5H11NO2/c1-6(2,3)4-5(7)8/h4H2,1-3H3"
        self.check(InchiLayers.INDEPENDENT_CHARGES, i1, i2)

    # D — Identity independent of double bond position
    def test_independent_double_bond_position(self):
        i1 = "InChI=1S/C18H34O2/c1-2-3-4-5-6-7-8-9-10-11-12-13-14-15-16-17-18(19)20/h9-10H,2-8,11-17H2,1H3,(H,19,20)/b10-9-"
        i2 = "InChI=1S/C18H34O2/c1-2-3-4-5-6-7-8-9-10-11-12-13-14-15-16-17-18(19)20/h8-9H,2-7,10-17H2,1H3,(H,19,20)/b9-8+"
        self.check(InchiLayers.INDEPENDENT_DOUBLE_BONDS, i1, i2)

    # E — Tautomeric identity
    def test_tautomeric_oxo_enol(self):
        i1 = "InChI=1S/C5H10O/c1-3-4-5(2)6/h3-4H2,1-2H3"
        i2 = "InChI=1S/C5H10O/c1-3-4-5(2)6/h4,6H,3H2,1-2H3/b5-4-"
        self.check(InchiLayers.TAUTOMERIC, i1, i2)

    def test_tautomeric_imine_enamine(self):
        i1 = "InChI=1S/C5H11N/c1-3-4-5(2)6/h6H,3-4H2,1-2H3"
        i2 = "InChI=1S/C5H11N/c1-3-4-5(2)6/h4H,3,6H2,1-2H3/b5-4-"
        self.check(InchiLayers.TAUTOMERIC, i1, i2)

    # F — Substituent position
    def test_stereochemical_identity_extended_smiles(self):
        i1 = "InChI=1S/C26H46NO8P/c1-6-8-10-12-13-14-15-17-19-26(29)35-24(22-32-25(28)18-16-11-9-7-2)23-34-36(30,31)33-21-20-27(3,4)5/h6-9,11,16,24H,10,12-15,17-23H2,1-5H3/b8-6-,9-7+,16-11+/t24-/m1/s1"
        i2 = "InChI=1S/C26H46NO8P/c1-5-6-7-8-9-10-11-16-19-26(30)33-22-24(17-14-12-13-15-18-25(28)29)23-35-36(31,32)34-21-20-27(2,3)4/h5-6,12-15,24H,7-11,16-23H2,1-4H3,(H-,28,29,31,32)/b6-5-,14-12+,15-13+/t24-/m1/s1"
        self.check(InchiLayers.STEREOCHEMICAL, i1, i2)

    # G — Identity independent of isotopes
    def test_independent_isotopes_hydrogen(self):
        i1 = "InChI=1S/C9H11NO2/c10-8(9(11)12)6-7-4-2-1-3-5-7/h1-5,8H,6,10H2,(H,11,12)/t8-/m0/s1"
        i2 = "InChI=1S/C9H11NO2/c10-8(9(11)12)6-7-4-2-1-3-5-7/h1-5,8H,6,10H2,(H,11,12)/t8-/m0/s1/i1D"
        self.check(InchiLayers.ISOTOPIC, i1, i2)

    def test_independent_isotopes_carbon(self):
        i1 = "InChI=1S/C9H11NO2/c10-8(9(11)12)6-7-4-2-1-3-5-7/h1-5,8H,6,10H2,(H,11,12)/t8-/m0/s1"
        i2 = "InChI=1S/C9H11NO2/c10-8(9(11)12)6-7-4-2-1-3-5-7/h1-5,8H,6,10H2,(H,11,12)/t8-/m0/s1/i9+1"
        self.check(InchiLayers.ISOTOPIC, i1, i2)
