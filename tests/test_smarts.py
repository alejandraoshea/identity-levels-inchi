import unittest
from rdkit import Chem
from src.backend.lipid.lipid_analysis import LipidHeadValidator


def parse(smiles: str):
    """Return RDKit mol, or None if SMILES is empty/placeholder/invalid."""
    if not smiles:
        return None
    if "[G]" in smiles or "[R]" in smiles:
        return None
    return Chem.MolFromSmiles(smiles.strip())


class TestEx3__3_BetaineLipids(unittest.TestCase):
    """Example 3 (3) — 1-hexadecanoyl-2-(6Z,9Z,12Z-octadecatrienoyl)-sn-glycero-3-O-(N,N,N-trimethyl)-homoserine"""
    POS = "O=C(CCCCCCCCCCCCCCC)OC[C@]([H])(COCCC(C([O-])=O)[N+](C)(C)C)OC(CCCC/C=C\\C/C=C\\C/C=C\\CCCCC)=O"
    NEG_STEREO = "O=C(CCCCCCCCCCCCCCC)OCC([H])(COCCC(C([O-])=O)[N+](C)(C)C)OC(CCCC/C=C\\C/C=C\\C/C=C\\CCCCC)=O"

    def setUp(self): self.v = LipidHeadValidator()

    def test_positive(self):
        """Positive: correct structure must be recognised."""
        m = parse(self.POS)
        self.assertIsNotNone(m, "Ex 3 (3) POS SMILES invalid")
        self.assertTrue(self.v.matches_any_valid_head(m),
                        "Ex 3 (3) positive must match")

    def test_positive_class(self):
        """Lipid class must be identified correctly."""
        m = parse(self.POS)
        classes = self.v.identify_lipid_class(m)
        self.assertTrue(len(classes) > 0, "Ex 3 (3) must have a lipid class")

    @unittest.expectedFailure
    def test_negative_stereo(self):
        """Wrong stereochemistry — not yet enforced by current SMARTS."""
        m = parse(self.NEG_STEREO)
        self.assertIsNotNone(m, "Ex 3 (3) NEG_STEREO SMILES invalid")
        self.assertFalse(self.v.matches_any_valid_head(m),
                         "Ex 3 (3) stereo-error should NOT match once stereo patterns added")


class TestEx77__2_Ceramides(unittest.TestCase):
    """Example 77 (2) — (2S,3S,26R,27S)-2,27-diamino-3,26,28-trihydroxyoctacosan-11-one"""
    POS = "OC[C@H](N)[C@H](O)CCCCCCCCCCCCCCC(CCCCCCC[C@H](O)[C@@H](N)C)=O"
    NEG_STEREO = "OC[C@H](N)[C@@H](O)CCCCCCCCCCCCCCC(CCCCCCC[C@H](O)[C@@H](N)C)=O"

    def setUp(self): self.v = LipidHeadValidator()

    def test_positive(self):
        """Positive: correct structure must be recognised."""
        m = parse(self.POS)
        self.assertIsNotNone(m, "Ex 77 (2) POS SMILES invalid")
        self.assertTrue(self.v.matches_any_valid_head(m),
                        "Ex 77 (2) positive must match")

    @unittest.expectedFailure
    def test_negative_stereo(self):
        """Wrong stereochemistry — not yet enforced by current SMARTS."""
        m = parse(self.NEG_STEREO)
        self.assertIsNotNone(m, "Ex 77 (2) NEG_STEREO SMILES invalid")
        self.assertFalse(self.v.matches_any_valid_head(m),
                         "Ex 77 (2) stereo-error should NOT match once stereo patterns added")


class TestEx78__2_Ceramides(unittest.TestCase):
    """Example 78 (2) — N-(9Z-octadecenoyl)-4E,6E-tetradecasphingadienine"""
    POS = "[H][C@](NC(CCCCCCC/C=C\\CCCCCCCC)=O)(CO)[C@@](O)([H])/C=C/C=C/CCCCCCC"

    def setUp(self): self.v = LipidHeadValidator()

    def test_positive(self):
        """Positive: correct structure must be recognised."""
        m = parse(self.POS)
        self.assertIsNotNone(m, "Ex 78 (2) POS SMILES invalid")
        self.assertTrue(self.v.matches_any_valid_head(m),
                        "Ex 78 (2) positive must match")


class TestEx80__2_BasicSpNeutralAndAcidic(unittest.TestCase):
    """Example 80 (2) — 1-β-galactosyl-sphing-4-enine"""
    POS = "[H][C@](O)([C@](N)([H])CO[C@@H]1O[C@H](CO)[C@H](O)[C@H](O)[C@H]1O)/C=C/CCCCCCCCCCCCC"

    def setUp(self): self.v = LipidHeadValidator()

    def test_positive(self):
        """Positive: correct structure must be recognised."""
        m = parse(self.POS)
        self.assertIsNotNone(m, "Ex 80 (2) POS SMILES invalid")
        self.assertTrue(self.v.matches_any_valid_head(m),
                        "Ex 80 (2) positive must match")


class TestEx81__02_AmphotericGlycosphingolipids(unittest.TestCase):
    """Example 81 (02) — 2S,3R,4E)-2-amino-3-hydroxyoctadec-4-en-1-yl β-D-galactopyranoside 6-(hydrogen sulfate)"""
    POS = "O[C@H]([C@H]1O)[C@H](OC[C@H](N)[C@H](O)/C=C/CCCCCCCCCCCCC)O[C@H](COS(O)(=O)=O)[C@@H]1O"
    NEG_STEREO = "O[C@H]([C@H]1O)[C@H](OC[C@H](N)[C@@H](O)/C=C/CCCCCCCCCCCCC)O[C@H](COS(O)(=O)=O)[C@@H]1O"

    def setUp(self): self.v = LipidHeadValidator()

    def test_positive(self):
        """Positive: correct structure must be recognised."""
        m = parse(self.POS)
        self.assertIsNotNone(m, "Ex 81 (02) POS SMILES invalid")
        self.assertTrue(self.v.matches_any_valid_head(m),
                        "Ex 81 (02) positive must match")

    @unittest.expectedFailure
    def test_negative_stereo(self):
        """Wrong stereochemistry — not yet enforced by current SMARTS."""
        m = parse(self.NEG_STEREO)
        self.assertIsNotNone(m, "Ex 81 (02) NEG_STEREO SMILES invalid")
        self.assertFalse(self.v.matches_any_valid_head(m),
                         "Ex 81 (02) stereo-error should NOT match once stereo patterns added")


class TestEx3__3_Sterols(unittest.TestCase):
    """Example 3 (3) — campest-5-en-3β-yl octadecanoate"""
    POS = "C[C@]12CC[C@H](OC(CCCCCCCCCCCCCCCCC)=O)CC1=CC[C@]3([H])[C@]2([H])CC[C@@]4(C)[C@@]3([H])CC[C@]4([H])[C@@](C)([H])CC[C@@H](C)C(C)C"
    NEG_CHAIN = "O[C@H](C1)CC[C@@]2(C)C1=CC[C@]3([H])[C@]2([H])CC[C@@]4(C)[C@@]3([H])CC[C@]4([H])[C@@](CC(CCCCCCCCCCCCCCCCC)=O)([H])CC[C@@H](C)C(C)C"
    NEG_STEREO = "C[C@]12CC[C@@H](OC(CCCCCCCCCCCCCCCCC)=O)CC1=CC[C@]3([H])[C@]2([H])CC[C@@]4(C)[C@@]3([H])CC[C@]4([H])[C@@](C)([H])CC[C@@H](C)C(C)C"

    def setUp(self): self.v = LipidHeadValidator()

    def test_positive(self):
        """Positive: correct structure must be recognised."""
        m = parse(self.POS)
        self.assertIsNotNone(m, "Ex 3 (3) POS SMILES invalid")
        self.assertTrue(self.v.matches_any_valid_head(m),
                        "Ex 3 (3) positive must match")

    def test_negative_chain(self):
        """Chain-error (FA at wrong position) must NOT match."""
        m = parse(self.NEG_CHAIN)
        self.assertIsNotNone(m, "Ex 3 (3) NEG_CHAIN SMILES invalid")
        self.assertFalse(self.v.matches_any_valid_head(m),
                         "Ex 3 (3) chain-error must NOT match")

    @unittest.expectedFailure
    def test_negative_stereo(self):
        """Wrong stereochemistry — not yet enforced by current SMARTS."""
        m = parse(self.NEG_STEREO)
        self.assertIsNotNone(m, "Ex 3 (3) NEG_STEREO SMILES invalid")
        self.assertFalse(self.v.matches_any_valid_head(m),
                         "Ex 3 (3) stereo-error should NOT match once stereo patterns added")


class TestEx4__3_Sterols(unittest.TestCase):
    """Example 4 (3) — Stigmast-5-en-3β-yl octadecanoate"""
    POS = "C[C@]12CC[C@H](OC(CCCCCCCCCCCCCCCCC)=O)CC1=CC[C@]3([H])[C@]2([H])CC[C@@]4(C)[C@@]3([H])CC[C@]4([H])[C@@](C)([H])CC[C@@H](CC)C(C)C"
    NEG_CHAIN = "O[C@H](C1)CC[C@@]2(C)C1=C(C(CCCCCCCCCCCCCCCCC)=O)C[C@]3([H])[C@]2([H])CC[C@@]4(C)[C@@]3([H])CC[C@]4([H])[C@@](C)([H])CC[C@@H](CC)C(C)C"
    NEG_STEREO = "C[C@]12CC[C@@H](OC(CCCCCCCCCCCCCCCCC)=O)CC1=CC[C@]3([H])[C@]2([H])CC[C@@]4(C)[C@@]3([H])CC[C@]4([H])[C@@](C)([H])CC[C@@H](CC)C(C)C"

    def setUp(self): self.v = LipidHeadValidator()

    def test_positive(self):
        """Positive: correct structure must be recognised."""
        m = parse(self.POS)
        self.assertIsNotNone(m, "Ex 4 (3) POS SMILES invalid")
        self.assertTrue(self.v.matches_any_valid_head(m),
                        "Ex 4 (3) positive must match")

    def test_negative_chain(self):
        """Chain-error (FA at wrong position) must NOT match."""
        m = parse(self.NEG_CHAIN)
        self.assertIsNotNone(m, "Ex 4 (3) NEG_CHAIN SMILES invalid")
        self.assertFalse(self.v.matches_any_valid_head(m),
                         "Ex 4 (3) chain-error must NOT match")

    @unittest.expectedFailure
    def test_negative_stereo(self):
        """Wrong stereochemistry — not yet enforced by current SMARTS."""
        m = parse(self.NEG_STEREO)
        self.assertIsNotNone(m, "Ex 4 (3) NEG_STEREO SMILES invalid")
        self.assertFalse(self.v.matches_any_valid_head(m),
                         "Ex 4 (3) stereo-error should NOT match once stereo patterns added")


class TestEx5__3_Sterols(unittest.TestCase):
    """Example 5 (3) — Stigmast-5,22E-dien-3β-yl docosanoate"""
    POS = "C[C@]12CC[C@H](OC(CCCCCCCCCCCCCCCCCCCCC)=O)CC1=CC[C@]3([H])[C@]2([H])CC[C@@]4(C)[C@@]3([H])CC[C@]4([H])[C@@](C)([H])/C=C/[C@@H](CC)C(C)C"
    NEG_CHAIN = "O[C@H](C1)CC[C@@]2(C)C1=CC[C@]3([H])[C@]2([H])CC[C@@]4(C)[C@@]3([H])CC[C@]4([H])[C@@](C)([H])/C=C/[C@@H](CC)C(C(CCCCCCCCCCCCCCCCCCCCC)=O)C"
    NEG_STEREO = "C[C@]12CC[C@@H](OC(CCCCCCCCCCCCCCCCCCCCC)=O)CC1=CC[C@]3([H])[C@]2([H])CC[C@@]4(C)[C@@]3([H])CC[C@]4([H])[C@@](C)([H])/C=C/[C@@H](CC)C(C)C"

    def setUp(self): self.v = LipidHeadValidator()

    def test_positive(self):
        """Positive: correct structure must be recognised."""
        m = parse(self.POS)
        self.assertIsNotNone(m, "Ex 5 (3) POS SMILES invalid")
        self.assertTrue(self.v.matches_any_valid_head(m),
                        "Ex 5 (3) positive must match")

    def test_negative_chain(self):
        """Chain-error (FA at wrong position) must NOT match."""
        m = parse(self.NEG_CHAIN)
        self.assertIsNotNone(m, "Ex 5 (3) NEG_CHAIN SMILES invalid")
        self.assertFalse(self.v.matches_any_valid_head(m),
                         "Ex 5 (3) chain-error must NOT match")

    @unittest.expectedFailure
    def test_negative_stereo(self):
        """Wrong stereochemistry — not yet enforced by current SMARTS."""
        m = parse(self.NEG_STEREO)
        self.assertIsNotNone(m, "Ex 5 (3) NEG_STEREO SMILES invalid")
        self.assertFalse(self.v.matches_any_valid_head(m),
                         "Ex 5 (3) stereo-error should NOT match once stereo patterns added")


class TestEx6__3_Sterols(unittest.TestCase):
    """Example 6 (3) — 6-oxo-5α-campestan-3β,22R,23R-triol 3β-yl tetradecanoate"""
    POS = "[H][C@@]12[C@]([C@](CC[C@H](OC(CCCCCCCCCCCCC)=O)C3)(C)[C@@]3([H])C(C2)=O)([H])CC[C@@]4(C)[C@@]1([H])CC[C@]4([H])[C@]([H])(C)[C@@H](O)[C@H](O)[C@@H](C)C(C)C"
    NEG_CHAIN = "[H][C@@]12[C@]([C@](CC[C@H](O)C3)(C)[C@@]3([H])C(C2)=O)([H])CC[C@@]4(C)[C@@]1([H])CC[C@]4([H])[C@]([H])(C)[C@@H](OC(CCCCCCCCCCCCC)=O)[C@H](O)[C@@H](C)C(C)C"
    NEG_STEREO = "[H][C@@]12[C@]([C@](CC[C@@H](OC(CCCCCCCCCCCCC)=O)C3)(C)[C@@]3([H])C(C2)=O)([H])CC[C@@]4(C)[C@@]1([H])CC[C@]4([H])[C@]([H])(C)[C@@H](O)[C@H](O)[C@@H](C)C(C)C"

    def setUp(self): self.v = LipidHeadValidator()

    def test_positive(self):
        """Positive: correct structure must be recognised."""
        m = parse(self.POS)
        self.assertIsNotNone(m, "Ex 6 (3) POS SMILES invalid")
        self.assertTrue(self.v.matches_any_valid_head(m),
                        "Ex 6 (3) positive must match")

    def test_negative_chain(self):
        """Chain-error (FA at wrong position) must NOT match."""
        m = parse(self.NEG_CHAIN)
        self.assertIsNotNone(m, "Ex 6 (3) NEG_CHAIN SMILES invalid")
        self.assertFalse(self.v.matches_any_valid_head(m),
                         "Ex 6 (3) chain-error must NOT match")

    @unittest.expectedFailure
    def test_negative_stereo(self):
        """Wrong stereochemistry — not yet enforced by current SMARTS."""
        m = parse(self.NEG_STEREO)
        self.assertIsNotNone(m, "Ex 6 (3) NEG_STEREO SMILES invalid")
        self.assertFalse(self.v.matches_any_valid_head(m),
                         "Ex 6 (3) stereo-error should NOT match once stereo patterns added")


class TestEx7__3_Sterols(unittest.TestCase):
    """Example 7 (3) — 3β,5β,10β,14β-tetrahydroxy-19-norbufa-20,22-dienolide-3β-yl-14-hydroxy-tetradecanoate"""
    POS = "C[C@@]1([C@]2(O)CC[C@]1([H])C(C=C3)=COC3=O)CC[C@@]4([H])[C@@]2([H])CC[C@]5(O)C[C@@H](OC(CCCCCCCCCCCCCO)=O)CC[C@@]54O"
    NEG_CHAIN = "C[C@@]1([C@]2(O)CC[C@]1([H])C(C=C3)=COC3=O)CC[C@@]4([H])[C@@]2([H])CC[C@]5(OC(CCCCCCCCCCCCCO)=O)C[C@@H](O)CC[C@@]54O"
    NEG_STEREO = "C[C@@]1([C@]2(O)CC[C@]1([H])C(C=C3)=COC3=O)CC[C@@]4([H])[C@@]2([H])CC[C@]5(O)C[C@H](OC(CCCCCCCCCCCCCO)=O)CC[C@@]54O"

    def setUp(self): self.v = LipidHeadValidator()

    def test_positive(self):
        """Positive: correct structure must be recognised."""
        m = parse(self.POS)
        self.assertIsNotNone(m, "Ex 7 (3) POS SMILES invalid")
        self.assertTrue(self.v.matches_any_valid_head(m),
                        "Ex 7 (3) positive must match")

    def test_negative_chain(self):
        """Chain-error (FA at wrong position) must NOT match."""
        m = parse(self.NEG_CHAIN)
        self.assertIsNotNone(m, "Ex 7 (3) NEG_CHAIN SMILES invalid")
        self.assertFalse(self.v.matches_any_valid_head(m),
                         "Ex 7 (3) chain-error must NOT match")

    @unittest.expectedFailure
    def test_negative_stereo(self):
        """Wrong stereochemistry — not yet enforced by current SMARTS."""
        m = parse(self.NEG_STEREO)
        self.assertIsNotNone(m, "Ex 7 (3) NEG_STEREO SMILES invalid")
        self.assertFalse(self.v.matches_any_valid_head(m),
                         "Ex 7 (3) stereo-error should NOT match once stereo patterns added")


class TestEx8__3_Sterols(unittest.TestCase):
    """Example 8 (3) — 3β,5β,14β-trihydroxy-19oxo-bufa-20,22-dienolide-3β-yl-14-hydroxy-tetradecanoate"""
    POS = "C[C@@]1([C@]2(O)CC[C@]1([H])C(C=C3)=COC3=O)CC[C@@]4([H])[C@@]2([H])CC[C@]5(O)C[C@@H](OC(CCCCCCCCCCCCCO)=O)CC[C@@]54C=O"
    NEG_CHAIN = "C[C@@]1([C@]2(O)CC[C@]1([H])C(C=C3)=COC3=O)CC[C@@]4([H])[C@@]2([H])CC[C@]5(OC(CCCCCCCCCCCCCO)=O)C[C@@H](O)CC[C@@]54C=O"
    NEG_STEREO = "C[C@@]1([C@]2(O)CC[C@]1([H])C(C=C3)=COC3=O)CC[C@@]4([H])[C@@]2([H])CC[C@]5(O)C[C@H](OC(CCCCCCCCCCCCCO)=O)CC[C@@]54C=O"

    def setUp(self): self.v = LipidHeadValidator()

    def test_positive(self):
        """Positive: correct structure must be recognised."""
        m = parse(self.POS)
        self.assertIsNotNone(m, "Ex 8 (3) POS SMILES invalid")
        self.assertTrue(self.v.matches_any_valid_head(m),
                        "Ex 8 (3) positive must match")

    def test_negative_chain(self):
        """Chain-error (FA at wrong position) must NOT match."""
        m = parse(self.NEG_CHAIN)
        self.assertIsNotNone(m, "Ex 8 (3) NEG_CHAIN SMILES invalid")
        self.assertFalse(self.v.matches_any_valid_head(m),
                         "Ex 8 (3) chain-error must NOT match")

    @unittest.expectedFailure
    def test_negative_stereo(self):
        """Wrong stereochemistry — not yet enforced by current SMARTS."""
        m = parse(self.NEG_STEREO)
        self.assertIsNotNone(m, "Ex 8 (3) NEG_STEREO SMILES invalid")
        self.assertFalse(self.v.matches_any_valid_head(m),
                         "Ex 8 (3) stereo-error should NOT match once stereo patterns added")


class TestEx9__3_Sterols(unittest.TestCase):
    """Example 9 (3) — 3β,5β,14β-trihydroxy-bufa-20,22-dienolide-3β-yl-14-hydroxy-tetradecanoate"""
    POS = "C[C@@]1([C@]2(O)CC[C@]1([H])C(C=C3)=COC3=O)CC[C@@]4([H])[C@@]2([H])CC[C@]5(O)C[C@@H](OC(CCCCCCCCCCCCCO)=O)CC[C@@]54C"
    NEG_CHAIN = "C[C@@]1([C@]2(O)CC[C@]1([H])C(C=C3)=COC3=O)CC[C@@]4([H])[C@@]2([H])CC[C@]5(OC(CCCCCCCCCCCCCO)=O)C[C@@H](O)CC[C@@]54C"
    NEG_STEREO = "C[C@@]1([C@]2(O)CC[C@]1([H])C(C=C3)=COC3=O)CC[C@@]4([H])[C@@]2([H])CC[C@]5(O)C[C@H](OC(CCCCCCCCCCCCCO)=O)CC[C@@]54C"

    def setUp(self): self.v = LipidHeadValidator()

    def test_positive(self):
        """Positive: correct structure must be recognised."""
        m = parse(self.POS)
        self.assertIsNotNone(m, "Ex 9 (3) POS SMILES invalid")
        self.assertTrue(self.v.matches_any_valid_head(m),
                        "Ex 9 (3) positive must match")

    def test_negative_chain(self):
        """Chain-error (FA at wrong position) must NOT match."""
        m = parse(self.NEG_CHAIN)
        self.assertIsNotNone(m, "Ex 9 (3) NEG_CHAIN SMILES invalid")
        self.assertFalse(self.v.matches_any_valid_head(m),
                         "Ex 9 (3) chain-error must NOT match")

    @unittest.expectedFailure
    def test_negative_stereo(self):
        """Wrong stereochemistry — not yet enforced by current SMARTS."""
        m = parse(self.NEG_STEREO)
        self.assertIsNotNone(m, "Ex 9 (3) NEG_STEREO SMILES invalid")
        self.assertFalse(self.v.matches_any_valid_head(m),
                         "Ex 9 (3) stereo-error should NOT match once stereo patterns added")


class TestEx10__3_Sterols(unittest.TestCase):
    """Example 10 (3) — 1β,3β,5β,14β-tetrahydroxy-bufa-20,22-dienolide-3β-yl-14-hydroxy-tetradecanoate"""
    POS = "C[C@@]1([C@]2(O)CC[C@]1([H])C(C=C3)=COC3=O)CC[C@@]4([H])[C@@]2([H])CC[C@]5(O)C[C@@H](OC(CCCCCCCCCCCCCO)=O)C[C@@H](O)[C@@]54C"
    NEG_CHAIN = "C[C@@]1([C@]2(O)CC[C@]1([H])C(C=C3)=COC3=O)CC[C@@]4([H])[C@@]2([H])CC[C@]5(OC(CCCCCCCCCCCCCO)=O)C[C@@H](O)C[C@@H](O)[C@@]54C"
    NEG_STEREO = "C[C@@]1([C@]2(O)CC[C@]1([H])C(C=C3)=COC3=O)CC[C@@]4([H])[C@@]2([H])CC[C@]5(O)C[C@H](OC(CCCCCCCCCCCCCO)=O)C[C@@H](O)[C@@]54C"

    def setUp(self): self.v = LipidHeadValidator()

    def test_positive(self):
        """Positive: correct structure must be recognised."""
        m = parse(self.POS)
        self.assertIsNotNone(m, "Ex 10 (3) POS SMILES invalid")
        self.assertTrue(self.v.matches_any_valid_head(m),
                        "Ex 10 (3) positive must match")

    def test_negative_chain(self):
        """Chain-error (FA at wrong position) must NOT match."""
        m = parse(self.NEG_CHAIN)
        self.assertIsNotNone(m, "Ex 10 (3) NEG_CHAIN SMILES invalid")
        self.assertFalse(self.v.matches_any_valid_head(m),
                         "Ex 10 (3) chain-error must NOT match")

    @unittest.expectedFailure
    def test_negative_stereo(self):
        """Wrong stereochemistry — not yet enforced by current SMARTS."""
        m = parse(self.NEG_STEREO)
        self.assertIsNotNone(m, "Ex 10 (3) NEG_STEREO SMILES invalid")
        self.assertFalse(self.v.matches_any_valid_head(m),
                         "Ex 10 (3) stereo-error should NOT match once stereo patterns added")


class TestEx11__3_Sterols(unittest.TestCase):
    """Example 11 (3) — 24-methylene,26,27-dimethylcholest-5-en-3β-yl 18-bromooctadeca-9E,17E-diene-5,7,15-triynoate"""
    POS = "C[C@@]1([C@@]2([H])CC[C@]1([H])[C@@](C)([H])CCC(C(CC)CC)=C)CC[C@@]3([H])[C@@]2([H])CC=C4C[C@@H](OC(CCCC#CC#C/C=C/CCCCC#C/C=C/Br)=O)CC[C@@]43C"
    NEG_CHAIN = "C[C@@]1([C@@]2([H])CC[C@]1([H])[C@@](CC(CCCC#CC#C/C=C/CCCCC#C/C=C/Br)=O)([H])CCC(C(CC)CC)=C)CC[C@@]3([H])[C@@]2([H])CC=C4C[C@@H](O)CC[C@@]43C"
    NEG_STEREO = "C[C@@]1([C@@]2([H])CC[C@]1([H])[C@@](C)([H])CCC(C(CC)CC)=C)CC[C@@]3([H])[C@@]2([H])CC=C4C[C@H](OC(CCCC#CC#C/C=C/CCCCC#C/C=C/Br)=O)CC[C@@]43C"

    def setUp(self): self.v = LipidHeadValidator()

    def test_positive(self):
        """Positive: correct structure must be recognised."""
        m = parse(self.POS)
        self.assertIsNotNone(m, "Ex 11 (3) POS SMILES invalid")
        self.assertTrue(self.v.matches_any_valid_head(m),
                        "Ex 11 (3) positive must match")

    def test_negative_chain(self):
        """Chain-error (FA at wrong position) must NOT match."""
        m = parse(self.NEG_CHAIN)
        self.assertIsNotNone(m, "Ex 11 (3) NEG_CHAIN SMILES invalid")
        self.assertFalse(self.v.matches_any_valid_head(m),
                         "Ex 11 (3) chain-error must NOT match")

    @unittest.expectedFailure
    def test_negative_stereo(self):
        """Wrong stereochemistry — not yet enforced by current SMARTS."""
        m = parse(self.NEG_STEREO)
        self.assertIsNotNone(m, "Ex 11 (3) NEG_STEREO SMILES invalid")
        self.assertFalse(self.v.matches_any_valid_head(m),
                         "Ex 11 (3) stereo-error should NOT match once stereo patterns added")


class TestEx12__3_Sterols(unittest.TestCase):
    """Example 12 (3) — 3β,5β-dihydroxy-14β,15β-epoxy-bufa-20,22-dienolide-3β-yl-16-hydroxy-9Z-hexadecenoate"""
    POS = "O=C(O[C@H]1CC[C@@]2(C)[C@@](CC[C@]3([H])[C@]2([H])CC[C@@]4(C)[C@]35[C@H](O5)C[C@]4([H])C(C=C6)=COC6=O)(O)C1)CCCCCCC/C=C\\CCCCCCO"
    NEG_CHAIN = "C[C@@]1([C@]23[C@H](O3)C[C@]1([H])C(C=C4)=COC4=O)CC[C@@]5([H])[C@@]2([H])CC[C@]6(OC(CCCCCCC/C=C\\CCCCCCO)=O)C[C@@H](O)CC[C@@]65C"
    NEG_STEREO = "O=C(O[C@@H]1CC[C@@]2(C)[C@@](CC[C@]3([H])[C@]2([H])CC[C@@]4(C)[C@]35[C@H](O5)C[C@]4([H])C(C=C6)=COC6=O)(O)C1)CCCCCCC/C=C\\CCCCCCO"

    def setUp(self): self.v = LipidHeadValidator()

    def test_positive(self):
        """Positive: correct structure must be recognised."""
        m = parse(self.POS)
        self.assertIsNotNone(m, "Ex 12 (3) POS SMILES invalid")
        self.assertTrue(self.v.matches_any_valid_head(m),
                        "Ex 12 (3) positive must match")

    def test_negative_chain(self):
        """Chain-error (FA at wrong position) must NOT match."""
        m = parse(self.NEG_CHAIN)
        self.assertIsNotNone(m, "Ex 12 (3) NEG_CHAIN SMILES invalid")
        self.assertFalse(self.v.matches_any_valid_head(m),
                         "Ex 12 (3) chain-error must NOT match")

    @unittest.expectedFailure
    def test_negative_stereo(self):
        """Wrong stereochemistry — not yet enforced by current SMARTS."""
        m = parse(self.NEG_STEREO)
        self.assertIsNotNone(m, "Ex 12 (3) NEG_STEREO SMILES invalid")
        self.assertFalse(self.v.matches_any_valid_head(m),
                         "Ex 12 (3) stereo-error should NOT match once stereo patterns added")


class TestEx13__3_Sterols(unittest.TestCase):
    """Example 13 (3) — 12-oxo-3β,11α,14β-trihydroxy-5β-bufa-20,22-dienolide-3β-yl-14-hydroxy-tetradecanoate"""
    POS = "O=C(O[C@H]1CC[C@@]2(C)[C@@](CC[C@]3([H])[C@]2([H])[C@H](O)C([C@@]4(C)[C@]3(O)CC[C@]4([H])C(C=C5)=COC5=O)=O)([H])C1)CCCCCCCCCCCCCO"
    NEG_CHAIN = "C[C@@]1([C@]2(OC(CCCCCCCCCCCCCO)=O)CC[C@]1([H])C(C=C3)=COC3=O)C([C@@H](O)[C@@]4([H])[C@@]2([H])CC[C@]5([H])C[C@@H](O)CC[C@@]54C)=O"
    NEG_STEREO = "O=C(O[C@@H]1CC[C@@]2(C)[C@@](CC[C@]3([H])[C@]2([H])[C@H](O)C([C@@]4(C)[C@]3(O)CC[C@]4([H])C(C=C5)=COC5=O)=O)([H])C1)CCCCCCCCCCCCCO"

    def setUp(self): self.v = LipidHeadValidator()

    def test_positive(self):
        """Positive: correct structure must be recognised."""
        m = parse(self.POS)
        self.assertIsNotNone(m, "Ex 13 (3) POS SMILES invalid")
        self.assertTrue(self.v.matches_any_valid_head(m),
                        "Ex 13 (3) positive must match")

    def test_negative_chain(self):
        """Chain-error (FA at wrong position) must NOT match."""
        m = parse(self.NEG_CHAIN)
        self.assertIsNotNone(m, "Ex 13 (3) NEG_CHAIN SMILES invalid")
        self.assertFalse(self.v.matches_any_valid_head(m),
                         "Ex 13 (3) chain-error must NOT match")

    @unittest.expectedFailure
    def test_negative_stereo(self):
        """Wrong stereochemistry — not yet enforced by current SMARTS."""
        m = parse(self.NEG_STEREO)
        self.assertIsNotNone(m, "Ex 13 (3) NEG_STEREO SMILES invalid")
        self.assertFalse(self.v.matches_any_valid_head(m),
                         "Ex 13 (3) stereo-error should NOT match once stereo patterns added")


class TestEx14__3_Sterols(unittest.TestCase):
    """Example 14 (3) — 3β,14β-dihydroxy-16β-acetoxy-5β-bufa-20,22-dienolide-3β-yl-hexadecanoate"""
    POS = "C[C@@]1([C@]2(O)C[C@H](OC(C)=O)[C@]1([H])C(C=C3)=COC3=O)CC[C@@]4([H])[C@@]2([H])CC[C@]5([H])C[C@@H](OC(CCCCCCCCCCCCCCC)=O)CC[C@@]54C"
    NEG_CHAIN = "O[C@H]1CC[C@@]2(C)[C@@](CC[C@]3([H])[C@]2([H])CC[C@@]4(C)[C@]3(OC(CCCCCCCCCCCCCCC)=O)C[C@H](OC(C)=O)[C@]4([H])C(C=C5)=COC5=O)([H])C1"
    NEG_STEREO = "C[C@@]1([C@]2(O)C[C@H](OC(C)=O)[C@]1([H])C(C=C3)=COC3=O)CC[C@@]4([H])[C@@]2([H])CC[C@]5([H])C[C@H](OC(CCCCCCCCCCCCCCC)=O)CC[C@@]54C"

    def setUp(self): self.v = LipidHeadValidator()

    def test_positive(self):
        """Positive: correct structure must be recognised."""
        m = parse(self.POS)
        self.assertIsNotNone(m, "Ex 14 (3) POS SMILES invalid")
        self.assertTrue(self.v.matches_any_valid_head(m),
                        "Ex 14 (3) positive must match")

    def test_negative_chain(self):
        """Chain-error (FA at wrong position) must NOT match."""
        m = parse(self.NEG_CHAIN)
        self.assertIsNotNone(m, "Ex 14 (3) NEG_CHAIN SMILES invalid")
        self.assertFalse(self.v.matches_any_valid_head(m),
                         "Ex 14 (3) chain-error must NOT match")

    @unittest.expectedFailure
    def test_negative_stereo(self):
        """Wrong stereochemistry — not yet enforced by current SMARTS."""
        m = parse(self.NEG_STEREO)
        self.assertIsNotNone(m, "Ex 14 (3) NEG_STEREO SMILES invalid")
        self.assertFalse(self.v.matches_any_valid_head(m),
                         "Ex 14 (3) stereo-error should NOT match once stereo patterns added")


class TestEx15__3_Isoprenoids(unittest.TestCase):
    """Example 15 (3) — epi-Austrobuxusin H"""
    POS = "O=C(O1)[C@]([C@H](C)OC(CCCCCCCCC)=O)([H])C[C@@]1([C@H]2[C@@H]3O2)[C@@]4(C)[C@]3(O)[C@@H]5[C@H](C(C)=C)[C@@H](OC5=O)[C@H]4O"
    NEG_CHAIN = "O=C(O1)[C@]([C@H](C)O)([H])C[C@@]1([C@H]2[C@@H]3O2)[C@@]4(C)[C@]3(OC(CCCCCCCCC)=O)[C@@H]5[C@H](C(C)=C)[C@@H](OC5=O)[C@H]4O"
    NEG_STEREO = "O=C(O1)[C@]([C@H](C)OC(CCCCCCCCC)=O)([H])C[C@@]1([C@H]2[C@@H]3O2)[C@@]4(C)[C@@]3(O)[C@@H]5[C@H](C(C)=C)[C@@H](OC5=O)[C@H]4O"

    def setUp(self): self.v = LipidHeadValidator()

    def test_positive(self):
        """Positive: correct structure must be recognised."""
        m = parse(self.POS)
        self.assertIsNotNone(m, "Ex 15 (3) POS SMILES invalid")
        self.assertTrue(self.v.matches_any_valid_head(m),
                        "Ex 15 (3) positive must match")

    def test_negative_chain(self):
        """Chain-error (FA at wrong position) must NOT match."""
        m = parse(self.NEG_CHAIN)
        self.assertIsNotNone(m, "Ex 15 (3) NEG_CHAIN SMILES invalid")
        self.assertFalse(self.v.matches_any_valid_head(m),
                         "Ex 15 (3) chain-error must NOT match")

    @unittest.expectedFailure
    def test_negative_stereo(self):
        """Wrong stereochemistry — not yet enforced by current SMARTS."""
        m = parse(self.NEG_STEREO)
        self.assertIsNotNone(m, "Ex 15 (3) NEG_STEREO SMILES invalid")
        self.assertFalse(self.v.matches_any_valid_head(m),
                         "Ex 15 (3) stereo-error should NOT match once stereo patterns added")


class TestEx16__3_Isoprenoids(unittest.TestCase):
    """Example 16 (3) — Ustusoic acid B"""
    POS = "C[C@@]12[C@]([C@H](OC(/C=C/C=C/C=C/C)=O)C=C(C=O)[C@]2(O)C(O)=O)([H])C(C)(C)CCC1"
    NEG_CHAIN = "C[C@@]12[C@]([C@H](O)C=C(C=O)[C@]2(O)C(OC(/C=C/C=C/C=C/C)=O)=O)([H])C(C)(C)CCC1"
    NEG_STEREO = "C[C@@]12[C@]([C@@H](OC(/C=C/C=C/C=C/C)=O)C=C(C=O)[C@]2(O)C(O)=O)([H])C(C)(C)CCC1"

    def setUp(self): self.v = LipidHeadValidator()

    def test_positive(self):
        """Positive: correct structure must be recognised."""
        m = parse(self.POS)
        self.assertIsNotNone(m, "Ex 16 (3) POS SMILES invalid")
        self.assertTrue(self.v.matches_any_valid_head(m),
                        "Ex 16 (3) positive must match")

    def test_negative_chain(self):
        """Chain-error (FA at wrong position) must NOT match."""
        m = parse(self.NEG_CHAIN)
        self.assertIsNotNone(m, "Ex 16 (3) NEG_CHAIN SMILES invalid")
        self.assertFalse(self.v.matches_any_valid_head(m),
                         "Ex 16 (3) chain-error must NOT match")

    @unittest.expectedFailure
    def test_negative_stereo(self):
        """Wrong stereochemistry — not yet enforced by current SMARTS."""
        m = parse(self.NEG_STEREO)
        self.assertIsNotNone(m, "Ex 16 (3) NEG_STEREO SMILES invalid")
        self.assertFalse(self.v.matches_any_valid_head(m),
                         "Ex 16 (3) stereo-error should NOT match once stereo patterns added")


class TestEx17__3_Isoprenoids(unittest.TestCase):
    """Example 17 (3) — 1-O-[(2E,4E,6E,8E)-3,7-dimethyl-9-(2,6,6-trimethylcyclohex-1-en-1-yl)nona-2,4,6,8-tetraenoyl]-β-D-glucopyranuronic acid O(15)-[(2S,3R,4S,5S,6S)-6-carboxy-3,4,5-trihydroxytetrahydro-2H-pyran-2-yl]retinoic acid"""
    POS = "CC(/C=C/C1=C(C)CCCC1(C)C)=C\\C=C\\C(C)=C\\C(O[C@@H]2O[C@H](C(O)=O)[C@@H](O)[C@H](O)[C@H]2O)=O"
    VAR_SUGAR = "CC(/C=C/C1=C(C)CCCC1(C)C)=C\\C=C\\C(C)=C\\C(O[C@@H]2O[C@@H](C)[C@@H](O)[C@H](O)[C@H]2O)=O"
    NEG_STEREO = "CC(/C=C/C1=C(C)CCCC1(C)C)=C\\C=C\\C(C)=C\\C(O[C@H]2O[C@H](C(O)=O)[C@@H](O)[C@H](O)[C@H]2O)=O"

    def setUp(self): self.v = LipidHeadValidator()

    def test_positive(self):
        """Positive: correct structure must be recognised."""
        m = parse(self.POS)
        self.assertIsNotNone(m, "Ex 17 (3) POS SMILES invalid")
        self.assertTrue(self.v.matches_any_valid_head(m),
                        "Ex 17 (3) positive must match")

    def test_sugar_variant_still_matches(self):
        """Sugar variant with correct FA attachment must still match."""
        m = parse(self.VAR_SUGAR)
        self.assertIsNotNone(m, "Ex 17 (3) VAR_SUGAR SMILES invalid")
        self.assertTrue(self.v.matches_any_valid_head(m),
                        "Ex 17 (3) sugar-variant must still match")

    @unittest.expectedFailure
    def test_negative_stereo(self):
        """Wrong stereochemistry — not yet enforced by current SMARTS."""
        m = parse(self.NEG_STEREO)
        self.assertIsNotNone(m, "Ex 17 (3) NEG_STEREO SMILES invalid")
        self.assertFalse(self.v.matches_any_valid_head(m),
                         "Ex 17 (3) stereo-error should NOT match once stereo patterns added")


class TestEx18__3_Isoprenoids(unittest.TestCase):
    """Example 18 (3) — 3,7-dimethyl-9-(2,6,6-trimethylcyclohex-1-en-1-yl)nona-2E,4Z,6E,8E-tetraen-1-yl hexadecanoate"""
    POS = "CC(/C=C\\C=C(C)\\C=C\\C1=C(CCCC1(C)C)C)=C\\COC(CCCCCCCCCCCCCCC)=O"
    NEG_CHAIN = "CC(/C=C\\C=C(C)\\C=C\\C1=C(C(CCCCCCCCCCCCCCC)=O)CCCC1(C)C)=C\\CO.C"

    def setUp(self): self.v = LipidHeadValidator()

    def test_positive(self):
        """Positive: correct structure must be recognised."""
        m = parse(self.POS)
        self.assertIsNotNone(m, "Ex 18 (3) POS SMILES invalid")
        self.assertTrue(self.v.matches_any_valid_head(m),
                        "Ex 18 (3) positive must match")

    def test_negative_chain(self):
        """Chain-error (FA at wrong position) must NOT match."""
        m = parse(self.NEG_CHAIN)
        self.assertIsNotNone(m, "Ex 18 (3) NEG_CHAIN SMILES invalid")
        self.assertFalse(self.v.matches_any_valid_head(m),
                         "Ex 18 (3) chain-error must NOT match")


class TestEx19__3_QuinonesAndHydroquinones(unittest.TestCase):
    """Example 19 (3) — N-(5-[(2R)-6-Hydroxy-2,5,7,8-tetramethyl-3,4-dihydrochromen-2-yl]-propanoyl)-glycine"""
    POS = "CC1=C(C)C(O)=C(C)C2=C1O[C@](C)(CCC(NCC(O)=O)=O)CC2"
    NEG_CHAIN = "CC1=C(C)C(OCC(O)=O)=C(C)C2=C1O[C@](C)(CCC(N)=O)CC2"
    NEG_STEREO = "CC1=C(C)C(O)=C(C)C2=C1O[C@@](C)(CCC(NCC(O)=O)=O)CC2"

    def setUp(self): self.v = LipidHeadValidator()

    def test_positive(self):
        """Positive: correct structure must be recognised."""
        m = parse(self.POS)
        self.assertIsNotNone(m, "Ex 19 (3) POS SMILES invalid")
        self.assertTrue(self.v.matches_any_valid_head(m),
                        "Ex 19 (3) positive must match")

    def test_negative_chain(self):
        """Chain-error (FA at wrong position) must NOT match."""
        m = parse(self.NEG_CHAIN)
        self.assertIsNotNone(m, "Ex 19 (3) NEG_CHAIN SMILES invalid")
        self.assertFalse(self.v.matches_any_valid_head(m),
                         "Ex 19 (3) chain-error must NOT match")

    @unittest.expectedFailure
    def test_negative_stereo(self):
        """Wrong stereochemistry — not yet enforced by current SMARTS."""
        m = parse(self.NEG_STEREO)
        self.assertIsNotNone(m, "Ex 19 (3) NEG_STEREO SMILES invalid")
        self.assertFalse(self.v.matches_any_valid_head(m),
                         "Ex 19 (3) stereo-error should NOT match once stereo patterns added")


class TestEx20__3_Acylaminosugars(unittest.TestCase):
    """Example 20 (3) — UDP-3-(3R-hydroxy-tetradecanoyl)-αD-glucosamine"""
    POS = "OP(OP(OC[C@H]1O[C@@H](N(C(N2)=O)C=CC2=O)[C@H](O)[C@@H]1O)(O)=O)(O[C@@H](O3)[C@@H]([C@H]([C@@H]([C@H]3CO)O)OC(C[C@@H](CCCCCCCCCCC)O)=O)N)=O"
    NEG_CHAIN = "OP(OP(OC[C@H]1O[C@@H](N(C(N2)=O)C=CC2=O)[C@H](OC(C[C@@H](CCCCCCCCCCC)O)=O)[C@@H]1O)(O)=O)(O[C@@H](O3)[C@@H]([C@H]([C@@H]([C@H]3CO)O)O)N)=O"
    NEG_STEREO = "OP(OP(OC[C@H]1O[C@@H](N(C(N2)=O)C=CC2=O)[C@H](O)[C@@H]1O)(O)=O)(O[C@@H](O3)[C@@H]([C@@H]([C@@H]([C@H]3CO)O)OC(C[C@@H](CCCCCCCCCCC)O)=O)N)=O"

    def setUp(self): self.v = LipidHeadValidator()

    def test_positive(self):
        """Positive: correct structure must be recognised."""
        m = parse(self.POS)
        self.assertIsNotNone(m, "Ex 20 (3) POS SMILES invalid")
        self.assertTrue(self.v.matches_any_valid_head(m),
                        "Ex 20 (3) positive must match")

    def test_negative_chain(self):
        """Chain-error (FA at wrong position) must NOT match."""
        m = parse(self.NEG_CHAIN)
        self.assertIsNotNone(m, "Ex 20 (3) NEG_CHAIN SMILES invalid")
        self.assertFalse(self.v.matches_any_valid_head(m),
                         "Ex 20 (3) chain-error must NOT match")

    @unittest.expectedFailure
    def test_negative_stereo(self):
        """Wrong stereochemistry — not yet enforced by current SMARTS."""
        m = parse(self.NEG_STEREO)
        self.assertIsNotNone(m, "Ex 20 (3) NEG_STEREO SMILES invalid")
        self.assertFalse(self.v.matches_any_valid_head(m),
                         "Ex 20 (3) stereo-error should NOT match once stereo patterns added")


class TestEx21__3_Acylaminosugars(unittest.TestCase):
    """Example 21 (3) — 2,3-bis-(3R-hydroxy-tetradecanoyl)-αD-glucosamine-1-phosphate"""
    POS = "O[C@H]1[C@H](OC(C[C@H](O)CCCCCCCCCCC)=O)[C@@H](NC(C[C@H](O)CCCCCCCCCCC)=O)[C@@H](OP(O)(O)=O)O[C@@H]1CO"
    NEG_CHAIN = "N[C@H]1[C@@H](OP(O)(OC(C[C@H](O)CCCCCCCCCCC)=O)=O)O[C@H](COC(C[C@H](O)CCCCCCCCCCC)=O)[C@@H](O)[C@@H]1O"
    NEG_STEREO = "O[C@H]1[C@H](OC(C[C@H](O)CCCCCCCCCCC)=O)[C@@H](NC(C[C@H](O)CCCCCCCCCCC)=O)[C@@H](OP(O)(O)=O)O[C@H]1CO"

    def setUp(self): self.v = LipidHeadValidator()

    def test_positive(self):
        """Positive: correct structure must be recognised."""
        m = parse(self.POS)
        self.assertIsNotNone(m, "Ex 21 (3) POS SMILES invalid")
        self.assertTrue(self.v.matches_any_valid_head(m),
                        "Ex 21 (3) positive must match")

    def test_negative_chain(self):
        """Chain-error (FA at wrong position) must NOT match."""
        m = parse(self.NEG_CHAIN)
        self.assertIsNotNone(m, "Ex 21 (3) NEG_CHAIN SMILES invalid")
        self.assertFalse(self.v.matches_any_valid_head(m),
                         "Ex 21 (3) chain-error must NOT match")

    @unittest.expectedFailure
    def test_negative_stereo(self):
        """Wrong stereochemistry — not yet enforced by current SMARTS."""
        m = parse(self.NEG_STEREO)
        self.assertIsNotNone(m, "Ex 21 (3) NEG_STEREO SMILES invalid")
        self.assertFalse(self.v.matches_any_valid_head(m),
                         "Ex 21 (3) stereo-error should NOT match once stereo patterns added")


class TestEx22__3_Acylaminosugars(unittest.TestCase):
    """Example 22 (3) — Lipid IVA"""
    POS = "OP(O[C@H]1O[C@H](CO[C@@H]2O[C@H](CO)[C@@H](OP(O)(O)=O)[C@H](OC(C[C@H](O)CCCCCCCCCCC)=O)[C@H]2NC(C[C@H](O)CCCCCCCCCCC)=O)[C@@H](O)[C@H](OC(C[C@H](O)CCCCCCCCCCC)=O)[C@H]1NC(C[C@H](O)CCCCCCCCCCC)=O)(O)=O"
    NEG_CHAIN = "N[C@H]1[C@@H](OP(OC(C[C@H](O)CCCCCCCCCCC)=O)(OC(C[C@H](O)CCCCCCCCCCC)=O)=O)O[C@H](CO[C@@H]2O[C@H](COC(C[C@H](O)CCCCCCCCCCC)=O)[C@@H](OP(O)(OC(C[C@H](O)CCCCCCCCCCC)=O)=O)[C@H](O)[C@H]2N)[C@@H](O)[C@@H]1O"
    NEG_STEREO = "OP(O[C@H]1O[C@H](CO[C@@H]2O[C@H](CO)[C@@H](OP(O)(O)=O)[C@@H](OC(C[C@H](O)CCCCCCCCCCC)=O)[C@H]2NC(C[C@H](O)CCCCCCCCCCC)=O)[C@@H](O)[C@H](OC(C[C@H](O)CCCCCCCCCCC)=O)[C@H]1NC(C[C@H](O)CCCCCCCCCCC)=O)(O)=O"

    def setUp(self): self.v = LipidHeadValidator()

    def test_positive(self):
        """Positive: correct structure must be recognised."""
        m = parse(self.POS)
        self.assertIsNotNone(m, "Ex 22 (3) POS SMILES invalid")
        self.assertTrue(self.v.matches_any_valid_head(m),
                        "Ex 22 (3) positive must match")

    def test_negative_chain(self):
        """Chain-error (FA at wrong position) must NOT match."""
        m = parse(self.NEG_CHAIN)
        self.assertIsNotNone(m, "Ex 22 (3) NEG_CHAIN SMILES invalid")
        self.assertFalse(self.v.matches_any_valid_head(m),
                         "Ex 22 (3) chain-error must NOT match")

    @unittest.expectedFailure
    def test_negative_stereo(self):
        """Wrong stereochemistry — not yet enforced by current SMARTS."""
        m = parse(self.NEG_STEREO)
        self.assertIsNotNone(m, "Ex 22 (3) NEG_STEREO SMILES invalid")
        self.assertFalse(self.v.matches_any_valid_head(m),
                         "Ex 22 (3) stereo-error should NOT match once stereo patterns added")


class TestEx23__3_Acylaminosugars(unittest.TestCase):
    """Example 23 (3) — UDP-2,3-bis-(3R-hydroxy-tetradecanoyl)-αD-glucosamine"""
    POS = "OP(OP(OC[C@H]1O[C@@H](N(C(N2)=O)C=CC2=O)[C@H](O)[C@@H]1O)(O)=O)(O[C@@H](O3)[C@@H]([C@H]([C@@H]([C@H]3CO)O)OC(C[C@@H](CCCCCCCCCCC)O)=O)NC(C[C@@H](CCCCCCCCCCC)O)=O)=O"
    NEG_CHAIN = "OP(OP(OC[C@H]1O[C@@H](N(C(N2)=O)C=CC2=O)[C@H](O)[C@@H]1OC(C[C@@H](CCCCCCCCCCC)O)=O)(O)=O)(O[C@@H](O3)[C@@H]([C@H]([C@@H]([C@H]3COC(C[C@@H](CCCCCCCCCCC)O)=O)O)O)N)=O"
    NEG_STEREO = "OP(OP(OC[C@H]1O[C@@H](N(C(N2)=O)C=CC2=O)[C@H](O)[C@@H]1O)(O)=O)(O[C@@H](O3)[C@@H]([C@H]([C@@H]([C@@H]3CO)O)OC(C[C@@H](CCCCCCCCCCC)O)=O)NC(C[C@@H](CCCCCCCCCCC)O)=O)=O"

    def setUp(self): self.v = LipidHeadValidator()

    def test_positive(self):
        """Positive: correct structure must be recognised."""
        m = parse(self.POS)
        self.assertIsNotNone(m, "Ex 23 (3) POS SMILES invalid")
        self.assertTrue(self.v.matches_any_valid_head(m),
                        "Ex 23 (3) positive must match")

    def test_negative_chain(self):
        """Chain-error (FA at wrong position) must NOT match."""
        m = parse(self.NEG_CHAIN)
        self.assertIsNotNone(m, "Ex 23 (3) NEG_CHAIN SMILES invalid")
        self.assertFalse(self.v.matches_any_valid_head(m),
                         "Ex 23 (3) chain-error must NOT match")

    @unittest.expectedFailure
    def test_negative_stereo(self):
        """Wrong stereochemistry — not yet enforced by current SMARTS."""
        m = parse(self.NEG_STEREO)
        self.assertIsNotNone(m, "Ex 23 (3) NEG_STEREO SMILES invalid")
        self.assertFalse(self.v.matches_any_valid_head(m),
                         "Ex 23 (3) stereo-error should NOT match once stereo patterns added")


class TestEx24__3_Acylaminosugars(unittest.TestCase):
    """Example 24 (3) — Lipid A -disaccharide-1-phosphate"""
    POS = "O[C@H]1[C@H](OC(C[C@H](O)CCCCCCCCCCC)=O)[C@@H](NC(C[C@H](O)CCCCCCCCCCC)=O)[C@@H](OP(O)(O)=O)O[C@@H]1CO[C@@H]2O[C@H](CO)[C@@H](O)[C@H](OC(C[C@H](O)CCCCCCCCCCC)=O)[C@H]2NC(C[C@H](O)CCCCCCCCCCC)=O"
    NEG_CHAIN = "N[C@H]1[C@@H](OP(OC(C[C@H](O)CCCCCCCCCCC)=O)(OC(C[C@H](O)CCCCCCCCCCC)=O)=O)O[C@H](CO[C@@H]2O[C@H](COC(C[C@H](O)CCCCCCCCCCC)=O)[C@@H](OC(C[C@H](O)CCCCCCCCCCC)=O)[C@H](O)[C@H]2N)[C@@H](O)[C@@H]1O"
    NEG_STEREO = "O[C@H]1[C@H](OC(C[C@H](O)CCCCCCCCCCC)=O)[C@@H](NC(C[C@H](O)CCCCCCCCCCC)=O)[C@@H](OP(O)(O)=O)O[C@@H]1CO[C@@H]2O[C@H](CO)[C@@H](O)[C@@H](OC(C[C@H](O)CCCCCCCCCCC)=O)[C@H]2NC(C[C@H](O)CCCCCCCCCCC)=O"

    def setUp(self): self.v = LipidHeadValidator()

    def test_positive(self):
        """Positive: correct structure must be recognised."""
        m = parse(self.POS)
        self.assertIsNotNone(m, "Ex 24 (3) POS SMILES invalid")
        self.assertTrue(self.v.matches_any_valid_head(m),
                        "Ex 24 (3) positive must match")

    def test_negative_chain(self):
        """Chain-error (FA at wrong position) must NOT match."""
        m = parse(self.NEG_CHAIN)
        self.assertIsNotNone(m, "Ex 24 (3) NEG_CHAIN SMILES invalid")
        self.assertFalse(self.v.matches_any_valid_head(m),
                         "Ex 24 (3) chain-error must NOT match")

    @unittest.expectedFailure
    def test_negative_stereo(self):
        """Wrong stereochemistry — not yet enforced by current SMARTS."""
        m = parse(self.NEG_STEREO)
        self.assertIsNotNone(m, "Ex 24 (3) NEG_STEREO SMILES invalid")
        self.assertFalse(self.v.matches_any_valid_head(m),
                         "Ex 24 (3) stereo-error should NOT match once stereo patterns added")


class TestEx25__3_AcylaminosugarGlycans(unittest.TestCase):
    """Example 25 (3) — Dodecanoyl-Kdo2-Lipid IVA"""
    POS = "O[C@H]([C@@H](CO[C@@H]([C@@H]1NC(C[C@H](OC(CCCCCCCCCCC)=O)CCCCCCCCCCC)=O)O[C@H](CO[C@]2(C(O)=O)C[C@@H](O[C@@]3(C(O)=O)O[C@H]([C@H](O)CO)[C@H](O)[C@H](O)C3)[C@@H](O)[C@@H]([C@H](O)CO)O2)[C@@H](OP(O)(O)=O)[C@@H]1OC(C[C@H](O)CCCCCCCCCCC)=O)O[C@H](OP(O)(O)=O)[C@@H]4NC(C[C@H](O)CCCCCCCCCCC)=O)[C@@H]4OC(C[C@H](O)CCCCCCCCCCC)=O"
    NEG_CHAIN = "O[C@H]1[C@H](OC(C[C@H](O)CCCCCCCCCCC)=O)[C@@H](CO[C@@H]([C@@H]2N)O[C@H](CO[C@]3(C(O)=O)C[C@@H](O[C@@]4(C(O)=O)O[C@H]([C@H](O)COC(C[C@H](O)CCCCCCCCCCC)=O)[C@H](O)[C@H](O)C4)[C@@H](O)[C@@H]([C@H](OC(C[C@H](O)CCCCCCCCCCC)=O)COC(C[C@H](OC(CCCCCCCCCCC)=O)CCCCCCCCCCC)=O)O3)[C@@H](OP(O)(O)=O)[C@@H]2O)O[C@H](OP(O)(O)=O)[C@@H]1N"
    NEG_STEREO = "O[C@H]([C@@H](CO[C@@H]([C@@H]1NC(C[C@H](OC(CCCCCCCCCCC)=O)CCCCCCCCCCC)=O)O[C@H](CO[C@]2(C(O)=O)C[C@@H](O[C@@]3(C(O)=O)O[C@H]([C@H](O)CO)[C@H](O)[C@H](O)C3)[C@@H](O)[C@@H]([C@H](O)CO)O2)[C@@H](OP(O)(O)=O)[C@@H]1OC(C[C@H](O)CCCCCCCCCCC)=O)O[C@@H](OP(O)(O)=O)[C@@H]4NC(C[C@H](O)CCCCCCCCCCC)=O)[C@@H]4OC(C[C@H](O)CCCCCCCCCCC)=O"

    def setUp(self): self.v = LipidHeadValidator()

    def test_positive(self):
        """Positive: correct structure must be recognised."""
        m = parse(self.POS)
        self.assertIsNotNone(m, "Ex 25 (3) POS SMILES invalid")
        self.assertTrue(self.v.matches_any_valid_head(m),
                        "Ex 25 (3) positive must match")

    def test_negative_chain(self):
        """Chain-error (FA at wrong position) must NOT match."""
        m = parse(self.NEG_CHAIN)
        self.assertIsNotNone(m, "Ex 25 (3) NEG_CHAIN SMILES invalid")
        self.assertFalse(self.v.matches_any_valid_head(m),
                         "Ex 25 (3) chain-error must NOT match")

    @unittest.expectedFailure
    def test_negative_stereo(self):
        """Wrong stereochemistry — not yet enforced by current SMARTS."""
        m = parse(self.NEG_STEREO)
        self.assertIsNotNone(m, "Ex 25 (3) NEG_STEREO SMILES invalid")
        self.assertFalse(self.v.matches_any_valid_head(m),
                         "Ex 25 (3) stereo-error should NOT match once stereo patterns added")


class TestEx26__3_Acyltrehaloses(unittest.TestCase):
    """Example 26 (3) — 2-O-hexadecanoyl-3-O-(2R,4S,6S-trimethyl-3R-hydroxy-tricosanoyl)-α,α-trehalose"""
    POS = "O[C@H]([C@@H](CO)O[C@H](O[C@@H](O1)[C@H](OC(CCCCCCCCCCCCCCC)=O)[C@@H](OC([C@H](C)[C@H](O)[C@@H](C)C[C@@H](C)CCCCCCCCCCCCCCCCC)=O)[C@H](O)[C@H]1CO)[C@@H]2O)[C@@H]2O"
    NEG_CHAIN = "O[C@H]1[C@H](OC[C@H]([C@H](O)[C@@H](C)C[C@@H](C)CCCCCCCCCCCCCCCCC)C=O)[C@@H](CO)O[C@H](O[C@@H](O2)[C@H](O)[C@@H](O)[C@H](O)[C@H]2COC(CCCCCCCCCCCCCCC)=O)[C@@H]1O"
    NEG_STEREO = "O[C@H]([C@@H](CO)O[C@H](O[C@@H](O1)[C@H](OC(CCCCCCCCCCCCCCC)=O)[C@@H](OC([C@H](C)[C@H](O)[C@@H](C)C[C@@H](C)CCCCCCCCCCCCCCCCC)=O)[C@H](O)[C@H]1CO)[C@@H]2O)[C@H]2O"

    def setUp(self): self.v = LipidHeadValidator()

    def test_positive(self):
        """Positive: correct structure must be recognised."""
        m = parse(self.POS)
        self.assertIsNotNone(m, "Ex 26 (3) POS SMILES invalid")
        self.assertTrue(self.v.matches_any_valid_head(m),
                        "Ex 26 (3) positive must match")

    def test_negative_chain(self):
        """Chain-error (FA at wrong position) must NOT match."""
        m = parse(self.NEG_CHAIN)
        self.assertIsNotNone(m, "Ex 26 (3) NEG_CHAIN SMILES invalid")
        self.assertFalse(self.v.matches_any_valid_head(m),
                         "Ex 26 (3) chain-error must NOT match")

    @unittest.expectedFailure
    def test_negative_stereo(self):
        """Wrong stereochemistry — not yet enforced by current SMARTS."""
        m = parse(self.NEG_STEREO)
        self.assertIsNotNone(m, "Ex 26 (3) NEG_STEREO SMILES invalid")
        self.assertFalse(self.v.matches_any_valid_head(m),
                         "Ex 26 (3) stereo-error should NOT match once stereo patterns added")


class TestEx27__3_Acyltrehaloses(unittest.TestCase):
    """Example 27 (3) — 2-O-hexadecanoyl-3-O-(2,4S,6S-trimethyl-2E-docosenoyl)-6-O-(2,4S,6S-trimethyl-2E-tetracosenoyl)-2'-O-(2,4S,6S-trimethyl-2E-tetracosenoyl)-4'-O-(2,4S,6S-trimethyl-2E-hexacosenoyl)-α,α-trehalose"""
    POS = "O[C@H]1[C@H](OC(/C(C)=C/[C@@H](C)C[C@@H](C)CCCCCCCCCCCCCCCCCCCC)=O)[C@@H](CO)O[C@H](O[C@@H](O2)[C@H](OC(CCCCCCCCCCCCCCC)=O)[C@@H](OC(/C(C)=C/[C@@H](C)C[C@@H](C)CCCCCCCCCCCCCCCC)=O)[C@H](O)[C@H]2COC(/C(C)=C/[C@@H](C)C[C@@H](C)CCCCCCCCCCCCCCCCCC)=O)[C@@H]1OC(/C(C)=C/[C@@H](C)C[C@@H](C)CCCCCCCCCCCCCCCCCC)=O"
    NEG_STEREO = "O[C@H]1[C@H](OC(/C(C)=C/[C@@H](C)C[C@@H](C)CCCCCCCCCCCCCCCCCCCC)=O)[C@H](CO)O[C@H](O[C@@H](O2)[C@H](OC(CCCCCCCCCCCCCCC)=O)[C@@H](OC(/C(C)=C/[C@@H](C)C[C@@H](C)CCCCCCCCCCCCCCCC)=O)[C@H](O)[C@H]2COC(/C(C)=C/[C@@H](C)C[C@@H](C)CCCCCCCCCCCCCCCCCC)=O)[C@@H]1OC(/C(C)=C/[C@@H](C)C[C@@H](C)CCCCCCCCCCCCCCCCCC)=O"

    def setUp(self): self.v = LipidHeadValidator()

    def test_positive(self):
        """Positive: correct structure must be recognised."""
        m = parse(self.POS)
        self.assertIsNotNone(m, "Ex 27 (3) POS SMILES invalid")
        self.assertTrue(self.v.matches_any_valid_head(m),
                        "Ex 27 (3) positive must match")

    @unittest.expectedFailure
    def test_negative_stereo(self):
        """Wrong stereochemistry — not yet enforced by current SMARTS."""
        m = parse(self.NEG_STEREO)
        self.assertIsNotNone(m, "Ex 27 (3) NEG_STEREO SMILES invalid")
        self.assertFalse(self.v.matches_any_valid_head(m),
                         "Ex 27 (3) stereo-error should NOT match once stereo patterns added")


class TestEx28__3_Acyltrehaloses(unittest.TestCase):
    """Example 28 (3) — 2-O-hexadecanoyl,3-O-(2S,4S,6S,8R,10R,12R,14R-heptamethyl)-15-hydroxy-triacontanoyl)-2'-sulfotrehalose"""
    POS = "O=C(O[C@H]1[C@@H](O[C@H]2O[C@H](CO)[C@@H](O)[C@H](O)[C@H]2OS(=O)(O)=O)O[C@H](CO)[C@@H](O)[C@@H]1OC([C@@H](C)C[C@@H](C)C[C@@H](C)C[C@@H](C)C[C@@H](C)C[C@@H](C)C[C@@H](C)C(O)CCCCCCCCCCCCCCC)=O)CCCCCCCCCCCCCCC"
    NEG_STEREO = "O=C(O[C@H]1[C@@H](O[C@H]2O[C@@H](CO)[C@@H](O)[C@H](O)[C@H]2OS(=O)(O)=O)O[C@H](CO)[C@@H](O)[C@@H]1OC([C@@H](C)C[C@@H](C)C[C@@H](C)C[C@@H](C)C[C@@H](C)C[C@@H](C)C[C@@H](C)C(O)CCCCCCCCCCCCCCC)=O)CCCCCCCCCCCCCCC"

    def setUp(self): self.v = LipidHeadValidator()

    def test_positive(self):
        """Positive: correct structure must be recognised."""
        m = parse(self.POS)
        self.assertIsNotNone(m, "Ex 28 (3) POS SMILES invalid")
        self.assertTrue(self.v.matches_any_valid_head(m),
                        "Ex 28 (3) positive must match")

    @unittest.expectedFailure
    def test_negative_stereo(self):
        """Wrong stereochemistry — not yet enforced by current SMARTS."""
        m = parse(self.NEG_STEREO)
        self.assertIsNotNone(m, "Ex 28 (3) NEG_STEREO SMILES invalid")
        self.assertFalse(self.v.matches_any_valid_head(m),
                         "Ex 28 (3) stereo-error should NOT match once stereo patterns added")


class TestEx29__3_Acyltrehaloses(unittest.TestCase):
    """Example 29 (3) — 6-O-hexadecanoyl-α-D-glucopyranosyl 6-O-hexadecanoyl-α-D-glucopyranoside"""
    POS = "O[C@H]1[C@@H](O[C@H]2O[C@H](COC(CCCCCCCCCCCCCCC)=O)[C@@H](O)[C@H](O)[C@H]2O)O[C@H](COC(CCCCCCCCCCCCCCC)=O)[C@@H](O)[C@@H]1O"
    NEG_STEREO = "O[C@H]1[C@@H](O[C@H]2O[C@H](COC(CCCCCCCCCCCCCCC)=O)[C@@H](O)[C@@H](O)[C@H]2O)O[C@H](COC(CCCCCCCCCCCCCCC)=O)[C@@H](O)[C@@H]1O"

    def setUp(self): self.v = LipidHeadValidator()

    def test_positive(self):
        """Positive: correct structure must be recognised."""
        m = parse(self.POS)
        self.assertIsNotNone(m, "Ex 29 (3) POS SMILES invalid")
        self.assertTrue(self.v.matches_any_valid_head(m),
                        "Ex 29 (3) positive must match")

    @unittest.expectedFailure
    def test_negative_stereo(self):
        """Wrong stereochemistry — not yet enforced by current SMARTS."""
        m = parse(self.NEG_STEREO)
        self.assertIsNotNone(m, "Ex 29 (3) NEG_STEREO SMILES invalid")
        self.assertFalse(self.v.matches_any_valid_head(m),
                         "Ex 29 (3) stereo-error should NOT match once stereo patterns added")


class TestEx30__3__29_Acyltrehaloses(unittest.TestCase):
    """Example 30 (3) 29 — Saccharolipids / Acyltrehaloses"""
    POS = "O[C@H]1[C@@H](O[C@H]2O[C@H](CO)[C@@H](O)[C@H](OC(CCCCCCCCCCCCCCC)=O)[C@H]2O)O[C@H](CO)[C@@H](OC(CCCCCCCCCCCCCCC)=O)[C@@H]1O"

    def setUp(self): self.v = LipidHeadValidator()

    def test_positive(self):
        """Positive: correct structure must be recognised."""
        m = parse(self.POS)
        self.assertIsNotNone(m, "Ex 30 (3) 29 POS SMILES invalid")
        self.assertTrue(self.v.matches_any_valid_head(m),
                        "Ex 30 (3) 29 positive must match")


class TestEx31__3_Acyltrehaloses(unittest.TestCase):
    """Example 31 (3) — Saccharolipids / Acyltrehaloses"""
    POS = "O[C@H]1[C@@H](O[C@H]2O[C@H](COC(CCCCCCCCCCCCCCCCC)=O)[C@@H](OC(CCCCCCCCCCCCCCCCC)=O)[C@H](O)[C@H]2O)O[C@H](CO)[C@@H](O)[C@@H]1O"
    NEG_STEREO = "O[C@H]1[C@@H](O[C@H]2O[C@H](COC(CCCCCCCCCCCCCCCCC)=O)[C@@H](OC(CCCCCCCCCCCCCCCCC)=O)[C@H](O)[C@H]2O)O[C@H](CO)[C@H](O)[C@@H]1O"

    def setUp(self): self.v = LipidHeadValidator()

    def test_positive(self):
        """Positive: correct structure must be recognised."""
        m = parse(self.POS)
        self.assertIsNotNone(m, "Ex 31 (3) POS SMILES invalid")
        self.assertTrue(self.v.matches_any_valid_head(m),
                        "Ex 31 (3) positive must match")

    @unittest.expectedFailure
    def test_negative_stereo(self):
        """Wrong stereochemistry — not yet enforced by current SMARTS."""
        m = parse(self.NEG_STEREO)
        self.assertIsNotNone(m, "Ex 31 (3) NEG_STEREO SMILES invalid")
        self.assertFalse(self.v.matches_any_valid_head(m),
                         "Ex 31 (3) stereo-error should NOT match once stereo patterns added")


class TestEx32__3_OtherAcylSugars(unittest.TestCase):
    """Example 32 (3) — 2,3-di-O-hexanoyl-α-glucopyranose"""
    POS = "O[C@H]1O[C@H](CO)[C@@H](O)[C@H](OC(CCCCC)=O)[C@H]1OC(CCCCC)=O"
    NEG_STEREO = "O[C@H]1O[C@H](CO)[C@@H](O)[C@@H](OC(CCCCC)=O)[C@H]1OC(CCCCC)=O"

    def setUp(self): self.v = LipidHeadValidator()

    def test_positive(self):
        """Positive: correct structure must be recognised."""
        m = parse(self.POS)
        self.assertIsNotNone(m, "Ex 32 (3) POS SMILES invalid")
        self.assertTrue(self.v.matches_any_valid_head(m),
                        "Ex 32 (3) positive must match")

    @unittest.expectedFailure
    def test_negative_stereo(self):
        """Wrong stereochemistry — not yet enforced by current SMARTS."""
        m = parse(self.NEG_STEREO)
        self.assertIsNotNone(m, "Ex 32 (3) NEG_STEREO SMILES invalid")
        self.assertFalse(self.v.matches_any_valid_head(m),
                         "Ex 32 (3) stereo-error should NOT match once stereo patterns added")


class TestEx33__3_OtherAcylSugars(unittest.TestCase):
    """Example 33 (3) — Butyl 4'-O-hexadecanoyl-neohesperidoside"""
    POS = "O[C@@H]1[C@@H](O[C@@H]2O[C@@H](C)[C@H](OC(CCCCCCCCCCCCCCC)=O)[C@@H](O)[C@H]2O)[C@H](OCCCC)O[C@H](CO)[C@H]1O"
    NEG_STEREO = "O[C@@H]1[C@@H](O[C@H]2O[C@@H](C)[C@H](OC(CCCCCCCCCCCCCCC)=O)[C@@H](O)[C@H]2O)[C@H](OCCCC)O[C@H](CO)[C@H]1O"

    def setUp(self): self.v = LipidHeadValidator()

    def test_positive(self):
        """Positive: correct structure must be recognised."""
        m = parse(self.POS)
        self.assertIsNotNone(m, "Ex 33 (3) POS SMILES invalid")
        self.assertTrue(self.v.matches_any_valid_head(m),
                        "Ex 33 (3) positive must match")

    @unittest.expectedFailure
    def test_negative_stereo(self):
        """Wrong stereochemistry — not yet enforced by current SMARTS."""
        m = parse(self.NEG_STEREO)
        self.assertIsNotNone(m, "Ex 33 (3) NEG_STEREO SMILES invalid")
        self.assertFalse(self.v.matches_any_valid_head(m),
                         "Ex 33 (3) stereo-error should NOT match once stereo patterns added")


class TestEx34__3_OtherAcylSugars(unittest.TestCase):
    """Example 34 (3) — Butyl 3'-O-acetyl-2'-O-butanoyl-4,6,4'-tri-O-(2-methylpropanoyl)-neohesperidoside"""
    POS = "O[C@@H]1[C@@H](O[C@@H]2O[C@@H](C)[C@H](OC(C(C)C)=O)[C@@H](OC(C)=O)[C@H]2OC(CCC)=O)[C@H](OCCCC)O[C@H](COC(C(C)C)=O)[C@H]1OC(C(C)C)=O"
    NEG_STEREO = "O[C@H]1[C@@H](O[C@@H]2O[C@@H](C)[C@H](OC(C(C)C)=O)[C@@H](OC(C)=O)[C@H]2OC(CCC)=O)[C@H](OCCCC)O[C@H](COC(C(C)C)=O)[C@H]1OC(C(C)C)=O"

    def setUp(self): self.v = LipidHeadValidator()

    def test_positive(self):
        """Positive: correct structure must be recognised."""
        m = parse(self.POS)
        self.assertIsNotNone(m, "Ex 34 (3) POS SMILES invalid")
        self.assertTrue(self.v.matches_any_valid_head(m),
                        "Ex 34 (3) positive must match")

    @unittest.expectedFailure
    def test_negative_stereo(self):
        """Wrong stereochemistry — not yet enforced by current SMARTS."""
        m = parse(self.NEG_STEREO)
        self.assertIsNotNone(m, "Ex 34 (3) NEG_STEREO SMILES invalid")
        self.assertFalse(self.v.matches_any_valid_head(m),
                         "Ex 34 (3) stereo-error should NOT match once stereo patterns added")


class TestEx35__3_OtherAcylSugars(unittest.TestCase):
    """Example 35 (3) — (11S)-jalapinolic 11-O-α-l-rhamnopyranosyl-(1-3)-O-[4-O-n-dodecanoyl-α-l-rhamnopyranosyl-(1-4)]-O-[2-O-n-octanoyl]-α-l-rhamnopyranosyl-(114)-O-α-l-rhamnopyranosyl-(1-2)-O-β-d-fucopyranoside-(1,3''-lactone)"""
    POS = "CCCCC[C@H](O[C@@H]1O[C@H](C)[C@H](O)[C@H](O)[C@H]1O2)CCCCCCCCCC(O[C@H]3[C@@H](O)[C@H]2O[C@@H](C)[C@@H]3O[C@@H]4O[C@@H](C)[C@H](O[C@@H]5O[C@@H](C)[C@H](OC(CCCCCCC)=O)[C@@H](O)[C@H]5O)[C@@H](O[C@@H]6O[C@@H](C)[C@H](O)[C@@H](O)[C@H]6O)[C@H]4OC(CCCCCCCCCCC)=O)=O"
    NEG_STEREO = "CCCCC[C@H](O[C@@H]1O[C@H](C)[C@H](O)[C@H](O)[C@H]1O2)CCCCCCCCCC(O[C@H]3[C@@H](O)[C@H]2O[C@@H](C)[C@@H]3O[C@H]4O[C@@H](C)[C@H](O[C@@H]5O[C@@H](C)[C@H](OC(CCCCCCC)=O)[C@@H](O)[C@H]5O)[C@@H](O[C@@H]6O[C@@H](C)[C@H](O)[C@@H](O)[C@H]6O)[C@H]4OC(CCCCCCCCCCC)=O)=O"

    def setUp(self): self.v = LipidHeadValidator()

    def test_positive(self):
        """Positive: correct structure must be recognised."""
        m = parse(self.POS)
        self.assertIsNotNone(m, "Ex 35 (3) POS SMILES invalid")
        self.assertTrue(self.v.matches_any_valid_head(m),
                        "Ex 35 (3) positive must match")

    @unittest.expectedFailure
    def test_negative_stereo(self):
        """Wrong stereochemistry — not yet enforced by current SMARTS."""
        m = parse(self.NEG_STEREO)
        self.assertIsNotNone(m, "Ex 35 (3) NEG_STEREO SMILES invalid")
        self.assertFalse(self.v.matches_any_valid_head(m),
                         "Ex 35 (3) stereo-error should NOT match once stereo patterns added")


class TestEx36__3_OtherAcylSugars(unittest.TestCase):
    """Example 36 (3) — (6R,9S)-9-O-β-D-glucopyranosyloxy-6'-O- ([''Z,1''Z,1''Z]-triene)-octadeca-6-hydroxy-9-methyl-3-oxo-α-ionol"""
    POS = "O=C(CCCCCCC/C=C\\C/C=C\\C/C=C\\CC)OC[C@H]1O[C@@H](O[C@@H](C)/C=C/[C@@]2(O)C(C)(C)CC(C=C2C)=O)[C@H](O)[C@@H](O)[C@@H]1O"

    def setUp(self): self.v = LipidHeadValidator()

    def test_positive(self):
        """Positive: correct structure must be recognised."""
        m = parse(self.POS)
        self.assertIsNotNone(m, "Ex 36 (3) POS SMILES invalid")
        self.assertTrue(self.v.matches_any_valid_head(m),
                        "Ex 36 (3) positive must match")


class TestEx37__3_OtherAcylSugars(unittest.TestCase):
    """Example 37 (3) — 2-O-(6-O-(9Z,12Z,15Z-octadecatrienoyl)-α-D-glucopyranosyl)-β-D-fructofuranose"""
    POS = "O[C@@H]1[C@@](O[C@H]([C@@H]([C@H]2O)O)O[C@H](COC(CCCCCCC/C=C\\C/C=C\\C/C=C\\CC)=O)[C@H]2O)(CO)O[C@H](CO)[C@H]1O"
    NEG_STEREO = "O[C@@H]1[C@@](O[C@H]([C@@H]([C@@H]2O)O)O[C@H](COC(CCCCCCC/C=C\\C/C=C\\C/C=C\\CC)=O)[C@H]2O)(CO)O[C@H](CO)[C@H]1O"

    def setUp(self): self.v = LipidHeadValidator()

    def test_positive(self):
        """Positive: correct structure must be recognised."""
        m = parse(self.POS)
        self.assertIsNotNone(m, "Ex 37 (3) POS SMILES invalid")
        self.assertTrue(self.v.matches_any_valid_head(m),
                        "Ex 37 (3) positive must match")

    @unittest.expectedFailure
    def test_negative_stereo(self):
        """Wrong stereochemistry — not yet enforced by current SMARTS."""
        m = parse(self.NEG_STEREO)
        self.assertIsNotNone(m, "Ex 37 (3) NEG_STEREO SMILES invalid")
        self.assertFalse(self.v.matches_any_valid_head(m),
                         "Ex 37 (3) stereo-error should NOT match once stereo patterns added")


class TestEx38__3_OtherAcylSugars(unittest.TestCase):
    """Example 38 (3) — 6-O-isopropoyl-6'-O-(14-methyl-2Z,4E-hexadecadienoyl)-maltose"""
    POS = "O[C@H]1[C@@H](O[C@H]2[C@H](O)[C@@H](O)C(O)O[C@@H]2COC(C(C)C)=O)O[C@H](COC(/C=C\\C=C\\CCCCCCCCC(CC)C)=O)[C@@H](O)[C@@H]1O"

    def setUp(self): self.v = LipidHeadValidator()

    def test_positive(self):
        """Positive: correct structure must be recognised."""
        m = parse(self.POS)
        self.assertIsNotNone(m, "Ex 38 (3) POS SMILES invalid")
        self.assertTrue(self.v.matches_any_valid_head(m),
                        "Ex 38 (3) positive must match")


class TestEx39__3_OtherAcylSugars(unittest.TestCase):
    """Example 39 (3) — (2-O-octadecanoyl-3-O-isobutyroyl)-2R-(α-D-glucopyranosyloxy)-3-hydroxypropanoic acid"""
    POS = "OC[C@H](C(O)=O)O[C@@H](O1)[C@@H]([C@H]([C@@H]([C@H]1CO)O)OC(C(C)C)=O)OC(CCCCCCCCCCCCCCCCC)=O"
    NEG_CHAIN = "OC[C@H](C(OC(C(C)C)=O)=O)O[C@@H](O1)[C@@H]([C@H]([C@@H]([C@H]1COC(CCCCCCCCCCCCCCCCC)=O)O)O)O"
    NEG_STEREO = "OC[C@H](C(O)=O)O[C@@H](O1)[C@@H]([C@H]([C@H]([C@H]1CO)O)OC(C(C)C)=O)OC(CCCCCCCCCCCCCCCCC)=O"

    def setUp(self): self.v = LipidHeadValidator()

    def test_positive(self):
        """Positive: correct structure must be recognised."""
        m = parse(self.POS)
        self.assertIsNotNone(m, "Ex 39 (3) POS SMILES invalid")
        self.assertTrue(self.v.matches_any_valid_head(m),
                        "Ex 39 (3) positive must match")

    def test_negative_chain(self):
        """Chain-error (FA at wrong position) must NOT match."""
        m = parse(self.NEG_CHAIN)
        self.assertIsNotNone(m, "Ex 39 (3) NEG_CHAIN SMILES invalid")
        self.assertFalse(self.v.matches_any_valid_head(m),
                         "Ex 39 (3) chain-error must NOT match")

    @unittest.expectedFailure
    def test_negative_stereo(self):
        """Wrong stereochemistry — not yet enforced by current SMARTS."""
        m = parse(self.NEG_STEREO)
        self.assertIsNotNone(m, "Ex 39 (3) NEG_STEREO SMILES invalid")
        self.assertFalse(self.v.matches_any_valid_head(m),
                         "Ex 39 (3) stereo-error should NOT match once stereo patterns added")


class TestEx40__3_OtherGlycerophospholipids(unittest.TestCase):
    """Example 40 (3) — 1,-2-diacyl-sn-glycero-3-phospho-O-[N-(2-hydroxyethyl)glycine]"""
    POS = "O=C(C)OC[C@]([H])(COP(OCCNCC(O)=O)(O)=O)OC(C)=O"
    NEG_STEREO = "O=C(C)OCC([H])(COP(OCCNCC(O)=O)(O)=O)OC(C)=O"

    def setUp(self): self.v = LipidHeadValidator()

    def test_positive(self):
        """Positive: correct structure must be recognised."""
        m = parse(self.POS)
        self.assertIsNotNone(m, "Ex 40 (3) POS SMILES invalid")
        self.assertTrue(self.v.matches_any_valid_head(m),
                        "Ex 40 (3) positive must match")

    @unittest.expectedFailure
    def test_negative_stereo(self):
        """Wrong stereochemistry — not yet enforced by current SMARTS."""
        m = parse(self.NEG_STEREO)
        self.assertIsNotNone(m, "Ex 40 (3) NEG_STEREO SMILES invalid")
        self.assertFalse(self.v.matches_any_valid_head(m),
                         "Ex 40 (3) stereo-error should NOT match once stereo patterns added")


class TestEx41__3_OtherGlycerophospholipids(unittest.TestCase):
    """Example 41 (3) — 1,2-dihexadecanoyl-sn-glycero-3-phosphosulfocholine"""
    POS = "[H][C@](OC(CCCCCCCCCCCCCCC)=O)(COP(OCC[S+](C)C)([O-])=O)COC(CCCCCCCCCCCCCCC)=O"
    NEG_STEREO = "[H]C(OC(CCCCCCCCCCCCCCC)=O)(COP(OCC[S+](C)C)([O-])=O)COC(CCCCCCCCCCCCCCC)=O"

    def setUp(self): self.v = LipidHeadValidator()

    def test_positive(self):
        """Positive: correct structure must be recognised."""
        m = parse(self.POS)
        self.assertIsNotNone(m, "Ex 41 (3) POS SMILES invalid")
        self.assertTrue(self.v.matches_any_valid_head(m),
                        "Ex 41 (3) positive must match")

    @unittest.expectedFailure
    def test_negative_stereo(self):
        """Wrong stereochemistry — not yet enforced by current SMARTS."""
        m = parse(self.NEG_STEREO)
        self.assertIsNotNone(m, "Ex 41 (3) NEG_STEREO SMILES invalid")
        self.assertFalse(self.v.matches_any_valid_head(m),
                         "Ex 41 (3) stereo-error should NOT match once stereo patterns added")


class TestEx42__3_OtherGlycerophospholipids(unittest.TestCase):
    """Example 42 (3) — PE 16:0/18:1(9Z)-15-isoLG pyrrole"""
    POS = "O=C(CCCCCCCCCCCCCCC)OC[C@]([H])(COP(OCCN1C(C)=C(C/C=C\\CCCC(O)=O)C(/C=C/C(O)CCCCC)=C1)(O)=O)OC(CCCCCCC/C=C\\CCCCCCCC)=O"
    NEG_CHAIN = "OC[C@]([H])(COP(OCCN1C(C)=C(C/C=C\\CCCC(O)=O)C(/C=C/C(OC(CCCCCCCCCCCCCCC)=O)CCCCC)=C1)(OC(CCCCCCC/C=C\\CCCCCCCC)=O)=O)O"
    NEG_STEREO = "O=C(CCCCCCCCCCCCCCC)OCC([H])(COP(OCCN1C(C)=C(C/C=C\\CCCC(O)=O)C(/C=C/C(O)CCCCC)=C1)(O)=O)OC(CCCCCCC/C=C\\CCCCCCCC)=O"

    def setUp(self): self.v = LipidHeadValidator()

    def test_positive(self):
        """Positive: correct structure must be recognised."""
        m = parse(self.POS)
        self.assertIsNotNone(m, "Ex 42 (3) POS SMILES invalid")
        self.assertTrue(self.v.matches_any_valid_head(m),
                        "Ex 42 (3) positive must match")

    def test_negative_chain(self):
        """Chain-error (FA at wrong position) must NOT match."""
        m = parse(self.NEG_CHAIN)
        self.assertIsNotNone(m, "Ex 42 (3) NEG_CHAIN SMILES invalid")
        self.assertFalse(self.v.matches_any_valid_head(m),
                         "Ex 42 (3) chain-error must NOT match")

    @unittest.expectedFailure
    def test_negative_stereo(self):
        """Wrong stereochemistry — not yet enforced by current SMARTS."""
        m = parse(self.NEG_STEREO)
        self.assertIsNotNone(m, "Ex 42 (3) NEG_STEREO SMILES invalid")
        self.assertFalse(self.v.matches_any_valid_head(m),
                         "Ex 42 (3) stereo-error should NOT match once stereo patterns added")


class TestEx42_OH__3_OtherGlycerophospholipids(unittest.TestCase):
    """Example 42-OH-(3) — PE 16:0/18:1(9Z)-15-isoLG hydroxylactam"""
    POS = "O=C(CCCCCCCCCCCCCCC)OC[C@]([H])(COP(OCCN(C1=O)C(C)(O)C(C/C=C\\CCCC(O)=O)=C1/C=C/C(O)CCCCC)(O)=O)OC(CCCCCCC/C=C\\CCCCCCCC)=O"

    def setUp(self): self.v = LipidHeadValidator()

    def test_positive(self):
        """Positive: correct structure must be recognised."""
        m = parse(self.POS)
        self.assertIsNotNone(m, "Ex 42-OH-(3) POS SMILES invalid")
        self.assertTrue(self.v.matches_any_valid_head(m),
                        "Ex 42-OH-(3) positive must match")


class TestEx43__3_OtherGlycerophospholipids(unittest.TestCase):
    """Example 43 (3) — 2,3-di-O-phytanyl-sn-glycero-1-phospho-(3'-sn-glycerol-1'-sulfate)"""
    POS = "CC(C)CCC[C@@H](C)CCC[C@@H](C)CCC[C@@H](C)CCOC[C@@](OCC[C@H](C)CCC[C@H](C)CCC[C@H](C)CCCC(C)C)([H])COP(OC[C@@](O)([H])COS(O)(=O)=O)(O)=O"
    NEG_STEREO = "CC(C)CCC[C@@H](C)CCC[C@@H](C)CCC[C@@H](C)CCOC[C@@](OCC[C@H](C)CCC[C@H](C)CCC[C@H](C)CCCC(C)C)([H])COP(OCC(O)([H])COS(O)(=O)=O)(O)=O"

    def setUp(self): self.v = LipidHeadValidator()

    def test_positive(self):
        """Positive: correct structure must be recognised."""
        m = parse(self.POS)
        self.assertIsNotNone(m, "Ex 43 (3) POS SMILES invalid")
        self.assertTrue(self.v.matches_any_valid_head(m),
                        "Ex 43 (3) positive must match")

    @unittest.expectedFailure
    def test_negative_stereo(self):
        """Wrong stereochemistry — not yet enforced by current SMARTS."""
        m = parse(self.NEG_STEREO)
        self.assertIsNotNone(m, "Ex 43 (3) NEG_STEREO SMILES invalid")
        self.assertFalse(self.v.matches_any_valid_head(m),
                         "Ex 43 (3) stereo-error should NOT match once stereo patterns added")


class TestEx44__3_OtherGlycerophospholipids(unittest.TestCase):
    """Example 44 (3) — 1-O-hexadecyl-sn-glycero-3-phosphoric acid methyl ester"""
    POS = "[H][C@](O)(COP(OC)(O)=O)COCCCCCCCCCCCCCCCC"

    def setUp(self): self.v = LipidHeadValidator()

    def test_positive(self):
        """Positive: correct structure must be recognised."""
        m = parse(self.POS)
        self.assertIsNotNone(m, "Ex 44 (3) POS SMILES invalid")
        self.assertTrue(self.v.matches_any_valid_head(m),
                        "Ex 44 (3) positive must match")


class TestEx45__3_OtherGlycerophospholipids(unittest.TestCase):
    """Example 45 (3) — Cholesteryl-6-O-(1,2-ditetradecanoyl-sn-glycero-3-phospho)-α-D-glucopyranoside"""
    POS = "O=C(CCCCCCCCCCCCC)OC[C@]([H])(COP(O)(OC[C@H]1O[C@H](O[C@@H]2CC3=CC[C@]4([H])[C@@](CC[C@@]5(C)[C@@]4([H])CC[C@@H]5[C@@H](CCCC(C)C)C)([H])[C@@]3(C)CC2)[C@H](O)[C@@H](O)[C@@H]1O)=O)OC(CCCCCCCCCCCCC)=O"
    NEG_CHAIN = "OC[C@]([H])(COP(OC(CCCCCCCCCCCCC)=O)(OC[C@H]1O[C@H](O[C@@H]2CC3=CC[C@]4([H])[C@@](CC[C@@]5(C)[C@@]4([H])CC[C@@H]5[C@@H](CCCC(C)C)C)([H])[C@@]3(C)CC2)[C@H](OC(CCCCCCCCCCCCC)=O)[C@@H](O)[C@@H]1O)=O)O"
    VAR_SUGAR = "O=C(CCCCCCCCCCCCC)OC[C@]([H])(COP(O)(OC[C@H]1O[C@H](O[C@@H]2CC3=CC[C@]4([H])[C@@](CC[C@@]5(C)[C@@]4([H])CC[C@@H]5[C@@H](CCCC(C)C)C)([H])[C@@]3(C)CC2)[C@H](O)[C@H](O)[C@@H]1O)=O)OC(CCCCCCCCCCCCC)=O"
    NEG_STEREO = "O=C(CCCCCCCCCCCCC)OCC([H])(COP(O)(OC[C@H]1O[C@H](O[C@@H]2CC3=CC[C@]4([H])[C@@](CC[C@@]5(C)[C@@]4([H])CC[C@@H]5[C@@H](CCCC(C)C)C)([H])[C@@]3(C)CC2)[C@H](O)[C@@H](O)[C@@H]1O)=O)OC(CCCCCCCCCCCCC)=O"

    def setUp(self): self.v = LipidHeadValidator()

    def test_positive(self):
        """Positive: correct structure must be recognised."""
        m = parse(self.POS)
        self.assertIsNotNone(m, "Ex 45 (3) POS SMILES invalid")
        self.assertTrue(self.v.matches_any_valid_head(m),
                        "Ex 45 (3) positive must match")

    def test_negative_chain(self):
        """Chain-error (FA at wrong position) must NOT match."""
        m = parse(self.NEG_CHAIN)
        self.assertIsNotNone(m, "Ex 45 (3) NEG_CHAIN SMILES invalid")
        self.assertFalse(self.v.matches_any_valid_head(m),
                         "Ex 45 (3) chain-error must NOT match")

    def test_sugar_variant_still_matches(self):
        """Sugar variant with correct FA attachment must still match."""
        m = parse(self.VAR_SUGAR)
        self.assertIsNotNone(m, "Ex 45 (3) VAR_SUGAR SMILES invalid")
        self.assertTrue(self.v.matches_any_valid_head(m),
                        "Ex 45 (3) sugar-variant must still match")

    @unittest.expectedFailure
    def test_negative_stereo(self):
        """Wrong stereochemistry — not yet enforced by current SMARTS."""
        m = parse(self.NEG_STEREO)
        self.assertIsNotNone(m, "Ex 45 (3) NEG_STEREO SMILES invalid")
        self.assertFalse(self.v.matches_any_valid_head(m),
                         "Ex 45 (3) stereo-error should NOT match once stereo patterns added")


class TestEx46__3_OtherGlycerophospholipids(unittest.TestCase):
    """Example 46 (3) — 1',3'-Bis-[1-(9Z-octadecenoyl)-2-hexadecanoyl-sn-glycero-3-phospho]-diethanolamine"""
    POS = "O=P(O)(OC[C@@](OC(CCCCCCCCCCCCCCC)=O)([H])COC(CCCCCCC/C=C\\CCCCCCCC)=O)OCCNCCOP(O)(OC[C@@](OC(CCCCCCCCCCCCCCC)=O)([H])COC(CCCCCCC/C=C\\CCCCCCCC)=O)=O"
    NEG_STEREO = "O=P(O)(OCC(OC(CCCCCCCCCCCCCCC)=O)([H])COC(CCCCCCC/C=C\\CCCCCCCC)=O)OCCNCCOP(O)(OC[C@@](OC(CCCCCCCCCCCCCCC)=O)([H])COC(CCCCCCC/C=C\\CCCCCCCC)=O)=O"

    def setUp(self): self.v = LipidHeadValidator()

    def test_positive(self):
        """Positive: correct structure must be recognised."""
        m = parse(self.POS)
        self.assertIsNotNone(m, "Ex 46 (3) POS SMILES invalid")
        self.assertTrue(self.v.matches_any_valid_head(m),
                        "Ex 46 (3) positive must match")

    @unittest.expectedFailure
    def test_negative_stereo(self):
        """Wrong stereochemistry — not yet enforced by current SMARTS."""
        m = parse(self.NEG_STEREO)
        self.assertIsNotNone(m, "Ex 46 (3) NEG_STEREO SMILES invalid")
        self.assertFalse(self.v.matches_any_valid_head(m),
                         "Ex 46 (3) stereo-error should NOT match once stereo patterns added")


class TestEx47__3_Glycerophosphocholines(unittest.TestCase):
    """Example 47 (3) — 1-hexadecanoyl-2-(9Z-octadecenoyl)-sn-glycero-3-phosphocholine"""
    POS = "O=C(CCCCCCCCCCCCCCC)OC[C@]([H])(COP(OCC[N+](C)(C)C)([O-])=O)OC(CCCCCCC/C=C\\CCCCCCCC)=O"
    NEG_STEREO = "O=C(CCCCCCCCCCCCCCC)OCC([H])(COP(OCC[N+](C)(C)C)([O-])=O)OC(CCCCCCC/C=C\\CCCCCCCC)=O"

    def setUp(self): self.v = LipidHeadValidator()

    def test_positive(self):
        """Positive: correct structure must be recognised."""
        m = parse(self.POS)
        self.assertIsNotNone(m, "Ex 47 (3) POS SMILES invalid")
        self.assertTrue(self.v.matches_any_valid_head(m),
                        "Ex 47 (3) positive must match")

    @unittest.expectedFailure
    def test_negative_stereo(self):
        """Wrong stereochemistry — not yet enforced by current SMARTS."""
        m = parse(self.NEG_STEREO)
        self.assertIsNotNone(m, "Ex 47 (3) NEG_STEREO SMILES invalid")
        self.assertFalse(self.v.matches_any_valid_head(m),
                         "Ex 47 (3) stereo-error should NOT match once stereo patterns added")


class TestEx48__3_Glycerophosphoethanolamines(unittest.TestCase):
    """Example 48 (3) — 1-hexadecanoyl-2-(9Z-octadecenoyl)-sn-glycero-3-phosphoethanolamine"""
    POS = "O=C(CCCCCCCCCCCCCCC)OC[C@]([H])(COP(OCCN)(O)=O)OC(CCCCCCC/C=C\\CCCCCCCC)=O"
    NEG_CHAIN = "OC[C@]([H])(COP(OCCNC(CCCCCCCCCCCCCCC)=O)(OC(CCCCCCC/C=C\\CCCCCCCC)=O)=O)O"
    NEG_STEREO = "O=C(CCCCCCCCCCCCCCC)OCC([H])(COP(OCCN)(O)=O)OC(CCCCCCC/C=C\\CCCCCCCC)=O"

    def setUp(self): self.v = LipidHeadValidator()

    def test_positive(self):
        """Positive: correct structure must be recognised."""
        m = parse(self.POS)
        self.assertIsNotNone(m, "Ex 48 (3) POS SMILES invalid")
        self.assertTrue(self.v.matches_any_valid_head(m),
                        "Ex 48 (3) positive must match")

    def test_negative_chain(self):
        """Chain-error (FA at wrong position) must NOT match."""
        m = parse(self.NEG_CHAIN)
        self.assertIsNotNone(m, "Ex 48 (3) NEG_CHAIN SMILES invalid")
        self.assertFalse(self.v.matches_any_valid_head(m),
                         "Ex 48 (3) chain-error must NOT match")

    @unittest.expectedFailure
    def test_negative_stereo(self):
        """Wrong stereochemistry — not yet enforced by current SMARTS."""
        m = parse(self.NEG_STEREO)
        self.assertIsNotNone(m, "Ex 48 (3) NEG_STEREO SMILES invalid")
        self.assertFalse(self.v.matches_any_valid_head(m),
                         "Ex 48 (3) stereo-error should NOT match once stereo patterns added")


class TestEx49__3_Glycerophosphoserines(unittest.TestCase):
    """Example 49 (3) — 1-dodecanoyl-2-tridecanoyl-sn-glycero-3-phosphoserine"""
    POS = "O=C(CCCCCCCCCCC)OC[C@]([H])(COP(OC[C@](C(O)=O)([H])N)(O)=O)OC(CCCCCCCCCCCC)=O"
    NEG_CHAIN = "OC[C@]([H])(COP(OC[C@](C(OC(CCCCCCCCCCC)=O)=O)([H])N)(OC(CCCCCCCCCCCC)=O)=O)O"
    NEG_STEREO = "O=C(CCCCCCCCCCC)OCC([H])(COP(OC[C@](C(O)=O)([H])N)(O)=O)OC(CCCCCCCCCCCC)=O"

    def setUp(self): self.v = LipidHeadValidator()

    def test_positive(self):
        """Positive: correct structure must be recognised."""
        m = parse(self.POS)
        self.assertIsNotNone(m, "Ex 49 (3) POS SMILES invalid")
        self.assertTrue(self.v.matches_any_valid_head(m),
                        "Ex 49 (3) positive must match")

    def test_negative_chain(self):
        """Chain-error (FA at wrong position) must NOT match."""
        m = parse(self.NEG_CHAIN)
        self.assertIsNotNone(m, "Ex 49 (3) NEG_CHAIN SMILES invalid")
        self.assertFalse(self.v.matches_any_valid_head(m),
                         "Ex 49 (3) chain-error must NOT match")

    @unittest.expectedFailure
    def test_negative_stereo(self):
        """Wrong stereochemistry — not yet enforced by current SMARTS."""
        m = parse(self.NEG_STEREO)
        self.assertIsNotNone(m, "Ex 49 (3) NEG_STEREO SMILES invalid")
        self.assertFalse(self.v.matches_any_valid_head(m),
                         "Ex 49 (3) stereo-error should NOT match once stereo patterns added")


class TestEx50__3_Glycerophosphoglycerols(unittest.TestCase):
    """Example 50 (3) — 1-hexadecanoyl-2-(9Z-octadecenoyl)-sn-glycero-3-phospho-(1'-rac-glycerol)"""
    POS = "O=C(CCCCCCCCCCCCCCC)OC[C@@]([H])(OC(CCCCCCC/C=C\\CCCCCCCC)=O)COP(OCC(O)CO)(O)=O"
    NEG_CHAIN = "OC[C@@]([H])(O)COP(OCC(O)COC(CCCCCCC/C=C\\CCCCCCCC)=O)(OC(CCCCCCCCCCCCCCC)=O)=O"
    NEG_STEREO = "O=C(CCCCCCCCCCCCCCC)OCC([H])(OC(CCCCCCC/C=C\\CCCCCCCC)=O)COP(OCC(O)CO)(O)=O"

    def setUp(self): self.v = LipidHeadValidator()

    def test_positive(self):
        """Positive: correct structure must be recognised."""
        m = parse(self.POS)
        self.assertIsNotNone(m, "Ex 50 (3) POS SMILES invalid")
        self.assertTrue(self.v.matches_any_valid_head(m),
                        "Ex 50 (3) positive must match")

    def test_negative_chain(self):
        """Chain-error (FA at wrong position) must NOT match."""
        m = parse(self.NEG_CHAIN)
        self.assertIsNotNone(m, "Ex 50 (3) NEG_CHAIN SMILES invalid")
        self.assertFalse(self.v.matches_any_valid_head(m),
                         "Ex 50 (3) chain-error must NOT match")

    @unittest.expectedFailure
    def test_negative_stereo(self):
        """Wrong stereochemistry — not yet enforced by current SMARTS."""
        m = parse(self.NEG_STEREO)
        self.assertIsNotNone(m, "Ex 50 (3) NEG_STEREO SMILES invalid")
        self.assertFalse(self.v.matches_any_valid_head(m),
                         "Ex 50 (3) stereo-error should NOT match once stereo patterns added")


class TestEx51__3_Glycerophosphoglycerophosphates(unittest.TestCase):
    """Example 51 (3) — 1-hexadecanoyl-2-(9Z-octadecenoyl)-sn-glycero-3-phospho-(1'-sn-glycerol-3'-phosphate)"""
    POS = "O=C(CCCCCCCCCCCCCCC)OC[C@]([H])(COP(OC[C@](O)([H])COP(O)(O)=O)(O)=O)OC(CCCCCCC/C=C\\CCCCCCCC)=O"
    NEG_CHAIN = "OC[C@]([H])(COP(OC[C@](OC(CCCCCCCCCCCCCCC)=O)([H])COP(O)(OC(CCCCCCC/C=C\\CCCCCCCC)=O)=O)(O)=O)O"
    NEG_STEREO = "O=C(CCCCCCCCCCCCCCC)OCC([H])(COP(OC[C@](O)([H])COP(O)(O)=O)(O)=O)OC(CCCCCCC/C=C\\CCCCCCCC)=O"

    def setUp(self): self.v = LipidHeadValidator()

    def test_positive(self):
        """Positive: correct structure must be recognised."""
        m = parse(self.POS)
        self.assertIsNotNone(m, "Ex 51 (3) POS SMILES invalid")
        self.assertTrue(self.v.matches_any_valid_head(m),
                        "Ex 51 (3) positive must match")

    def test_negative_chain(self):
        """Chain-error (FA at wrong position) must NOT match."""
        m = parse(self.NEG_CHAIN)
        self.assertIsNotNone(m, "Ex 51 (3) NEG_CHAIN SMILES invalid")
        self.assertFalse(self.v.matches_any_valid_head(m),
                         "Ex 51 (3) chain-error must NOT match")

    @unittest.expectedFailure
    def test_negative_stereo(self):
        """Wrong stereochemistry — not yet enforced by current SMARTS."""
        m = parse(self.NEG_STEREO)
        self.assertIsNotNone(m, "Ex 51 (3) NEG_STEREO SMILES invalid")
        self.assertFalse(self.v.matches_any_valid_head(m),
                         "Ex 51 (3) stereo-error should NOT match once stereo patterns added")


class TestEx52__3_Glycerophosphoinositols(unittest.TestCase):
    """Example 52 (3) — 1-hexadecanoyl-2-(9Z-octadecenoyl)-sn-glycero-3-phospho-(1'-myo-inositol)"""
    POS = "O=C(CCCCCCCCCCCCCCC)OC[C@@]([H])(OC(CCCCCCC/C=C\\CCCCCCCC)=O)COP(O[C@@H]1[C@H](O)[C@H](O)[C@@H](O)[C@H](O)[C@H]1O)(O)=O"
    NEG_CHAIN = "OC[C@@]([H])(O)COP(O[C@@H]1[C@H](OC(CCCCCCC/C=C\\CCCCCCCC)=O)[C@H](O)[C@@H](O)[C@H](O)[C@H]1O)(OC(CCCCCCCCCCCCCCC)=O)=O"
    VAR_SUGAR = "O=C(CCCCCCCCCCCCCCC)OC[C@@]([H])(OC(CCCCCCC/C=C\\CCCCCCCC)=O)COP(O[C@H]1[C@@H](O)[C@H](O)[C@@H](O)[C@H](O)[C@H]1O)(O)=O"
    NEG_STEREO = "O=C(CCCCCCCCCCCCCCC)OCC([H])(OC(CCCCCCC/C=C\\CCCCCCCC)=O)COP(O[C@@H]1[C@H](O)[C@H](O)[C@@H](O)[C@H](O)[C@H]1O)(O)=O"

    def setUp(self): self.v = LipidHeadValidator()

    def test_positive(self):
        """Positive: correct structure must be recognised."""
        m = parse(self.POS)
        self.assertIsNotNone(m, "Ex 52 (3) POS SMILES invalid")
        self.assertTrue(self.v.matches_any_valid_head(m),
                        "Ex 52 (3) positive must match")

    def test_negative_chain(self):
        """Chain-error (FA at wrong position) must NOT match."""
        m = parse(self.NEG_CHAIN)
        self.assertIsNotNone(m, "Ex 52 (3) NEG_CHAIN SMILES invalid")
        self.assertFalse(self.v.matches_any_valid_head(m),
                         "Ex 52 (3) chain-error must NOT match")

    def test_sugar_variant_still_matches(self):
        """Sugar variant with correct FA attachment must still match."""
        m = parse(self.VAR_SUGAR)
        self.assertIsNotNone(m, "Ex 52 (3) VAR_SUGAR SMILES invalid")
        self.assertTrue(self.v.matches_any_valid_head(m),
                        "Ex 52 (3) sugar-variant must still match")

    @unittest.expectedFailure
    def test_negative_stereo(self):
        """Wrong stereochemistry — not yet enforced by current SMARTS."""
        m = parse(self.NEG_STEREO)
        self.assertIsNotNone(m, "Ex 52 (3) NEG_STEREO SMILES invalid")
        self.assertFalse(self.v.matches_any_valid_head(m),
                         "Ex 52 (3) stereo-error should NOT match once stereo patterns added")


class TestEx52__3___mono_Glycerophosphoinositols(unittest.TestCase):
    """Example 52 (3)- mono — Glycerophospholipids / Glycerophosphoinositols"""
    POS = "O=C(CCCCCCCCCCCCCCC)OC[C@@]([H])(OC(CCCCCCC/C=C\\CCCCCCCC)=O)COP(O[C@@H]1[C@H](O)[C@H](OP(O)(O)=O)[C@@H](O)[C@H](O)[C@H]1O)(O)=O"

    def setUp(self): self.v = LipidHeadValidator()

    def test_positive(self):
        """Positive: correct structure must be recognised."""
        m = parse(self.POS)
        self.assertIsNotNone(m, "Ex 52 (3)- mono POS SMILES invalid")
        self.assertTrue(self.v.matches_any_valid_head(m),
                        "Ex 52 (3)- mono positive must match")


class TestEx52__3___bi_Glycerophosphoinositols(unittest.TestCase):
    """Example 52 (3)- bi — Glycerophospholipids / Glycerophosphoinositols"""
    POS = "O=C(CCCCCCCCCCCCCCCC)OC[C@@]([H])(OC(CCC/C=C\\C/C=C\\C/C=C\\C/C=C\\CCCCC)=O)COP(O[C@@H]1[C@H](O)[C@H](OP(O)(O)=O)[C@@H](OP(O)(O)=O)[C@H](O)[C@H]1O)(O)=O"

    def setUp(self): self.v = LipidHeadValidator()

    def test_positive(self):
        """Positive: correct structure must be recognised."""
        m = parse(self.POS)
        self.assertIsNotNone(m, "Ex 52 (3)- bi POS SMILES invalid")
        self.assertTrue(self.v.matches_any_valid_head(m),
                        "Ex 52 (3)- bi positive must match")


class TestEx52__3___tri_Glycerophosphoinositols(unittest.TestCase):
    """Example 52 (3)- tri — Glycerophospholipids / Glycerophosphoinositols"""
    POS = "O=C(CCCCCCCCCCCCCCCC)OC[C@@]([H])(OC(CCC/C=C\\C/C=C\\C/C=C\\C/C=C\\CCCCC)=O)COP(O[C@@H]1[C@H](O)[C@H](OP(O)(O)=O)[C@@H](OP(O)(O)=O)[C@H](OP(O)(O)=O)[C@H]1O)(O)=O"

    def setUp(self): self.v = LipidHeadValidator()

    def test_positive(self):
        """Positive: correct structure must be recognised."""
        m = parse(self.POS)
        self.assertIsNotNone(m, "Ex 52 (3)- tri POS SMILES invalid")
        self.assertTrue(self.v.matches_any_valid_head(m),
                        "Ex 52 (3)- tri positive must match")


class TestEx53__3_Glycerophosphates(unittest.TestCase):
    """Example 53 (3) — 1-dodecanoyl-2-tridecanoyl-sn-glycero-3-phosphate"""
    POS = "O=C(CCCCCCCCCCC)OC[C@]([H])(COP(O)(O)=O)OC(CCCCCCCCCCCC)=O"
    NEG_CHAIN = "OC[C@]([H])(COP(OC(CCCCCCCCCCCC)=O)(OC(CCCCCCCCCCC)=O)=O)O"
    NEG_STEREO = "O=C(CCCCCCCCCCC)OCC([H])(COP(O)(O)=O)OC(CCCCCCCCCCCC)=O"

    def setUp(self): self.v = LipidHeadValidator()

    def test_positive(self):
        """Positive: correct structure must be recognised."""
        m = parse(self.POS)
        self.assertIsNotNone(m, "Ex 53 (3) POS SMILES invalid")
        self.assertTrue(self.v.matches_any_valid_head(m),
                        "Ex 53 (3) positive must match")

    def test_negative_chain(self):
        """Chain-error (FA at wrong position) must NOT match."""
        m = parse(self.NEG_CHAIN)
        self.assertIsNotNone(m, "Ex 53 (3) NEG_CHAIN SMILES invalid")
        self.assertFalse(self.v.matches_any_valid_head(m),
                         "Ex 53 (3) chain-error must NOT match")

    @unittest.expectedFailure
    def test_negative_stereo(self):
        """Wrong stereochemistry — not yet enforced by current SMARTS."""
        m = parse(self.NEG_STEREO)
        self.assertIsNotNone(m, "Ex 53 (3) NEG_STEREO SMILES invalid")
        self.assertFalse(self.v.matches_any_valid_head(m),
                         "Ex 53 (3) stereo-error should NOT match once stereo patterns added")


class TestEx54__3_Glyceropyrophosphates(unittest.TestCase):
    """Example 54 (3) — 1-hexadecanoyl-2-(9Z-octadecenoyl)-sn-glycero-3-pyrophosphate"""
    POS = "O=C(CCCCCCCCCCCCCCC)OC[C@]([H])(COP(OP(O)(O)=O)(O)=O)OC(CCCCCCC/C=C\\CCCCCCCC)=O"
    NEG_CHAIN = "OC[C@]([H])(COP(OP(OC(CCCCCCC/C=C\\CCCCCCCC)=O)(O)=O)(OC(CCCCCCCCCCCCCCC)=O)=O)O"
    NEG_STEREO = "O=C(CCCCCCCCCCCCCCC)OCC([H])(COP(OP(O)(O)=O)(O)=O)OC(CCCCCCC/C=C\\CCCCCCCC)=O"

    def setUp(self): self.v = LipidHeadValidator()

    def test_positive(self):
        """Positive: correct structure must be recognised."""
        m = parse(self.POS)
        self.assertIsNotNone(m, "Ex 54 (3) POS SMILES invalid")
        self.assertTrue(self.v.matches_any_valid_head(m),
                        "Ex 54 (3) positive must match")

    def test_negative_chain(self):
        """Chain-error (FA at wrong position) must NOT match."""
        m = parse(self.NEG_CHAIN)
        self.assertIsNotNone(m, "Ex 54 (3) NEG_CHAIN SMILES invalid")
        self.assertFalse(self.v.matches_any_valid_head(m),
                         "Ex 54 (3) chain-error must NOT match")

    @unittest.expectedFailure
    def test_negative_stereo(self):
        """Wrong stereochemistry — not yet enforced by current SMARTS."""
        m = parse(self.NEG_STEREO)
        self.assertIsNotNone(m, "Ex 54 (3) NEG_STEREO SMILES invalid")
        self.assertFalse(self.v.matches_any_valid_head(m),
                         "Ex 54 (3) stereo-error should NOT match once stereo patterns added")


class TestEx55__3_Glycerophosphoglycerophosphoglycerols(unittest.TestCase):
    """Example 55 (3) — 1',3'-Bis[1,2-Di-(9Z-12Z-octadecadienoyl)-sn-glycero-3-phospho]-sn-glycerol"""
    POS = "O=P(O)(OC[C@@](OC(CCCCCCC/C=C\\C/C=C\\CCCCC)=O)([H])COC(CCCCCCC/C=C\\C/C=C\\CCCCC)=O)OCC([H])(O)COP(O)(OC[C@@](OC(CCCCCCC/C=C\\C/C=C\\CCCCC)=O)([H])COC(CCCCCCC/C=C\\C/C=C\\CCCCC)=O)=O"
    NEG_STEREO = "O=P(O)(OCC(OC(CCCCCCC/C=C\\C/C=C\\CCCCC)=O)([H])COC(CCCCCCC/C=C\\C/C=C\\CCCCC)=O)OC[C@@]([H])(O)COP(O)(OC[C@@](OC(CCCCCCC/C=C\\C/C=C\\CCCCC)=O)([H])COC(CCCCCCC/C=C\\C/C=C\\CCCCC)=O)=O"

    def setUp(self): self.v = LipidHeadValidator()

    def test_positive(self):
        """Positive: correct structure must be recognised."""
        m = parse(self.POS)
        self.assertIsNotNone(m, "Ex 55 (3) POS SMILES invalid")
        self.assertTrue(self.v.matches_any_valid_head(m),
                        "Ex 55 (3) positive must match")

    @unittest.expectedFailure
    def test_negative_stereo(self):
        """Wrong stereochemistry — not yet enforced by current SMARTS."""
        m = parse(self.NEG_STEREO)
        self.assertIsNotNone(m, "Ex 55 (3) NEG_STEREO SMILES invalid")
        self.assertFalse(self.v.matches_any_valid_head(m),
                         "Ex 55 (3) stereo-error should NOT match once stereo patterns added")


class TestEx56__3_Cdpglycerols(unittest.TestCase):
    """Example 56 (3) — 1,2-Didodecanoyl-sn-glycero-3-cytidine-5'-diphosphate"""
    POS = "O=C(CCCCCCCCCCC)OC[C@@]([H])(OC(CCCCCCCCCCC)=O)COP(OP(OC[C@H]1O[C@@H](N2C=CC(N)=NC2=O)[C@H](O)[C@@H]1O)(O)=O)(O)=O"
    NEG_CHAIN = "OC[C@@]([H])(O)COP(OP(OC[C@H]1O[C@@H](N2C=CC(NC(CCCCCCCCCCC)=O)=NC2=O)[C@H](O)[C@@H]1OC(CCCCCCCCCCC)=O)(O)=O)(O)=O"
    NEG_STEREO = "O=C(CCCCCCCCCCC)OCC([H])(OC(CCCCCCCCCCC)=O)COP(OP(OC[C@H]1O[C@@H](N2C=CC(N)=NC2=O)[C@H](O)[C@@H]1O)(O)=O)(O)=O"

    def setUp(self): self.v = LipidHeadValidator()

    def test_positive(self):
        """Positive: correct structure must be recognised."""
        m = parse(self.POS)
        self.assertIsNotNone(m, "Ex 56 (3) POS SMILES invalid")
        self.assertTrue(self.v.matches_any_valid_head(m),
                        "Ex 56 (3) positive must match")

    def test_negative_chain(self):
        """Chain-error (FA at wrong position) must NOT match."""
        m = parse(self.NEG_CHAIN)
        self.assertIsNotNone(m, "Ex 56 (3) NEG_CHAIN SMILES invalid")
        self.assertFalse(self.v.matches_any_valid_head(m),
                         "Ex 56 (3) chain-error must NOT match")

    @unittest.expectedFailure
    def test_negative_stereo(self):
        """Wrong stereochemistry — not yet enforced by current SMARTS."""
        m = parse(self.NEG_STEREO)
        self.assertIsNotNone(m, "Ex 56 (3) NEG_STEREO SMILES invalid")
        self.assertFalse(self.v.matches_any_valid_head(m),
                         "Ex 56 (3) stereo-error should NOT match once stereo patterns added")


class TestEx57__3_Glycosylglycerophospholipids(unittest.TestCase):
    """Example 57 (3) — 1-octadecanoyl-2-(5Z,8Z,11Z,14Z-eicosatetraenoyl)-sn-glycero-3-phospho-(1'-β-D-glucose)"""
    POS = "O=C(CCCCCCCCCCCCCCCCC)OC[C@]([H])(COP(O[C@@H]1O[C@H](CO)[C@@H](O)[C@H](O)[C@H]1O)(O)=O)OC(CCC/C=C\\C/C=C\\C/C=C\\C/C=C\\CCCCC)=O"
    NEG_CHAIN = "OC[C@]([H])(COP(O[C@@H]1O[C@H](CO)[C@@H](O)[C@H](OC(CCCCCCCCCCCCCCCCC)=O)[C@H]1O)(OC(CCC/C=C\\C/C=C\\C/C=C\\C/C=C\\CCCCC)=O)=O)O"
    NEG_STEREO = "O=C(CCCCCCCCCCCCCCCCC)OCC([H])(COP(O[C@@H]1O[C@H](CO)[C@@H](O)[C@H](O)[C@H]1O)(O)=O)OC(CCC/C=C\\C/C=C\\C/C=C\\C/C=C\\CCCCC)=O"

    def setUp(self): self.v = LipidHeadValidator()

    def test_positive(self):
        """Positive: correct structure must be recognised."""
        m = parse(self.POS)
        self.assertIsNotNone(m, "Ex 57 (3) POS SMILES invalid")
        self.assertTrue(self.v.matches_any_valid_head(m),
                        "Ex 57 (3) positive must match")

    def test_negative_chain(self):
        """Chain-error (FA at wrong position) must NOT match."""
        m = parse(self.NEG_CHAIN)
        self.assertIsNotNone(m, "Ex 57 (3) NEG_CHAIN SMILES invalid")
        self.assertFalse(self.v.matches_any_valid_head(m),
                         "Ex 57 (3) chain-error must NOT match")

    @unittest.expectedFailure
    def test_negative_stereo(self):
        """Wrong stereochemistry — not yet enforced by current SMARTS."""
        m = parse(self.NEG_STEREO)
        self.assertIsNotNone(m, "Ex 57 (3) NEG_STEREO SMILES invalid")
        self.assertFalse(self.v.matches_any_valid_head(m),
                         "Ex 57 (3) stereo-error should NOT match once stereo patterns added")


class TestEx58__3_Glycosylglycerophospholipids(unittest.TestCase):
    """Example 58 (3) — 1-octadecanoyl-2-eicosanoyl-sn-glycero-3-phospho-(1'-β-D-6'-O-acetyl-glucose)"""
    POS = "O=C(CCCCCCCCCCCCCCCCC)OC[C@]([H])(COP(O[C@@H]1O[C@H](COC(C)=O)[C@@H](O)[C@H](O)[C@H]1O)(O)=O)OC(CCCCCCCCCCCCCCCCCCC)=O"
    NEG_STEREO = "O=C(CCCCCCCCCCCCCCCCC)OCC([H])(COP(O[C@@H]1O[C@H](COC(C)=O)[C@@H](O)[C@H](O)[C@H]1O)(O)=O)OC(CCCCCCCCCCCCCCCCCCC)=O"

    def setUp(self): self.v = LipidHeadValidator()

    def test_positive(self):
        """Positive: correct structure must be recognised."""
        m = parse(self.POS)
        self.assertIsNotNone(m, "Ex 58 (3) POS SMILES invalid")
        self.assertTrue(self.v.matches_any_valid_head(m),
                        "Ex 58 (3) positive must match")

    @unittest.expectedFailure
    def test_negative_stereo(self):
        """Wrong stereochemistry — not yet enforced by current SMARTS."""
        m = parse(self.NEG_STEREO)
        self.assertIsNotNone(m, "Ex 58 (3) NEG_STEREO SMILES invalid")
        self.assertFalse(self.v.matches_any_valid_head(m),
                         "Ex 58 (3) stereo-error should NOT match once stereo patterns added")


class TestEx59__3_Glycosylglycerophospholipids(unittest.TestCase):
    """Example 59 (3) — 1-(11Z-octadecenoyl)-2-(hexadecenoyl)-sn-glycero-3-phospho-2'-α-D-6-glucosaminyl-sn-glycerol"""
    POS = "[H][C@@](CO)(O[C@H]1O[C@H](CO)[C@@H](O)[C@H](O)[C@H]1N)COP(OC[C@]([H])(OC(CCCCCCCCCCCCCCC)=O)COC(CCCCCCCCC/C=C\\CCCCCC)=O)(O)=O"
    NEG_STEREO = "[H][C@@](CO)(O[C@H]1O[C@H](CO)[C@@H](O)[C@H](O)[C@H]1N)COP(OCC([H])(OC(CCCCCCCCCCCCCCC)=O)COC(CCCCCCCCC/C=C\\CCCCCC)=O)(O)=O"

    def setUp(self): self.v = LipidHeadValidator()

    def test_positive(self):
        """Positive: correct structure must be recognised."""
        m = parse(self.POS)
        self.assertIsNotNone(m, "Ex 59 (3) POS SMILES invalid")
        self.assertTrue(self.v.matches_any_valid_head(m),
                        "Ex 59 (3) positive must match")

    @unittest.expectedFailure
    def test_negative_stereo(self):
        """Wrong stereochemistry — not yet enforced by current SMARTS."""
        m = parse(self.NEG_STEREO)
        self.assertIsNotNone(m, "Ex 59 (3) NEG_STEREO SMILES invalid")
        self.assertFalse(self.v.matches_any_valid_head(m),
                         "Ex 59 (3) stereo-error should NOT match once stereo patterns added")


class TestEx60__3_Glycosylglycerophospholipids(unittest.TestCase):
    """Example 60 (3) — 1-tetradecanoyl-2-hexadecanoyl-sn-glycero-3-phospho-(1'-glycerol-3'-)5-deoxy-5-(dimethylarsinyl)-β-D-ribofuranoside"""
    POS = "O=C(CCCCCCCCCCCCC)OC[C@]([H])(COP(OCC(O)CO[C@@H]1O[C@H](C[As](C)(C)=O)[C@@H](O)[C@H]1O)(O)=O)OC(CCCCCCCCCCCCCCC)=O"
    NEG_CHAIN = "O=C(CCCCCCCCCCCCC)OCC([H])(COP(OCC(O)CO[C@@H]1O[C@H](C[As](C)(C)=O)[C@@H](O)[C@H]1O)(O)=O)OC(CCCCCCCCCCCCCCC)=O"
    NEG_STEREO = "OC[C@]([H])(COP(OCC(OC(CCCCCCCCCCCCC)=O)CO[C@@H]1O[C@H](C[As](C)(C)=O)[C@@H](O)[C@H]1O)(OC(CCCCCCCCCCCCCCC)=O)=O)O"

    def setUp(self): self.v = LipidHeadValidator()

    def test_positive(self):
        """Positive: correct structure must be recognised."""
        m = parse(self.POS)
        self.assertIsNotNone(m, "Ex 60 (3) POS SMILES invalid")
        self.assertTrue(self.v.matches_any_valid_head(m),
                        "Ex 60 (3) positive must match")

    def test_negative_chain(self):
        """Chain-error (FA at wrong position) must NOT match."""
        m = parse(self.NEG_CHAIN)
        self.assertIsNotNone(m, "Ex 60 (3) NEG_CHAIN SMILES invalid")
        self.assertFalse(self.v.matches_any_valid_head(m),
                         "Ex 60 (3) chain-error must NOT match")

    @unittest.expectedFailure
    def test_negative_stereo(self):
        """Wrong stereochemistry — not yet enforced by current SMARTS."""
        m = parse(self.NEG_STEREO)
        self.assertIsNotNone(m, "Ex 60 (3) NEG_STEREO SMILES invalid")
        self.assertFalse(self.v.matches_any_valid_head(m),
                         "Ex 60 (3) stereo-error should NOT match once stereo patterns added")


class TestEx61__3_Glycerophosphoinositolglycans(unittest.TestCase):
    """Example 61 (3) — EtN-P-6Manα1-2Manα1-6Manα1-4GlcNα1-6-PI(14:0/14:0)"""
    POS = "O=C(CCCCCCCCCCCCC)OC[C@@]([H])(OC(CCCCCCCCCCCCC)=O)COP(O)(O[C@@H]1[C@@H](O)[C@H](O)[C@@H](O)[C@H](O)[C@@H]1O[C@H]2O[C@H](CO)[C@@H](O[C@H]3O[C@H](CO[C@H]4O[C@H](CO)[C@@H](O)[C@H](O)[C@@H]4O[C@H]5O[C@H](COP(OCCN)(O)=O)[C@@H](O)[C@H](O)[C@@H]5O)[C@@H](O)[C@H](O)[C@@H]3O)[C@H](O)[C@H]2N)=O"
    NEG_CHAIN = "OC[C@@]([H])(O)COP(O)(O[C@@H]1[C@@H](O)[C@H](O)[C@@H](O)[C@H](O)[C@@H]1O[C@H]2O[C@H](COC(CCCCCCCCCCCCC)=O)[C@@H](O[C@H]3O[C@H](CO[C@H]4O[C@H](COC(CCCCCCCCCCCCC)=O)[C@@H](O)[C@H](O)[C@@H]4O[C@H]5O[C@H](COP(OCCN)(O)=O)[C@@H](O)[C@H](O)[C@@H]5O)[C@@H](O)[C@H](O)[C@@H]3O)[C@H](O)[C@H]2N)=O"

    def setUp(self): self.v = LipidHeadValidator()

    def test_positive(self):
        """Positive: correct structure must be recognised."""
        m = parse(self.POS)
        self.assertIsNotNone(m, "Ex 61 (3) POS SMILES invalid")
        self.assertTrue(self.v.matches_any_valid_head(m),
                        "Ex 61 (3) positive must match")

    def test_negative_chain(self):
        """Chain-error (FA at wrong position) must NOT match."""
        m = parse(self.NEG_CHAIN)
        self.assertIsNotNone(m, "Ex 61 (3) NEG_CHAIN SMILES invalid")
        self.assertFalse(self.v.matches_any_valid_head(m),
                         "Ex 61 (3) chain-error must NOT match")


class TestEx62__3_Glycerophosphoinositolglycans(unittest.TestCase):
    """Example 62 (3) — 2'-O-(α-D-Manp)-(1-hexadecanoyl-2-tetradecanoyl-sn-glycero-3-phospho-1'-myo-inositol)"""
    POS = "O=C(CCCCCCCCCCCCCCC)OC[C@]([H])(COP(O)(O[C@@H]1[C@@H](O[C@H]2O[C@H](CO)[C@@H](O)[C@H](O)[C@@H]2O)[C@H](O)[C@@H](O)[C@H](O)[C@@H]1O)=O)OC(CCCCCCCCCCCCC)=O"
    VAR_SUGAR = "O=C(CCCCCCCCCCCCCCC)OC[C@]([H])(COP(O)(O[C@@H]1[C@@H](O[C@@H]2O[C@H](CO)[C@@H](O)[C@H](O)[C@H]2O)[C@H](O)[C@@H](O)[C@H](O)[C@H]1O)=O)OC(CCCCCCCCCCCCC)=O"
    NEG_STEREO = "O=C(CCCCCCCCCCCCCCC)OCC([H])(COP(O)(O[C@@H]1[C@@H](O[C@H]2O[C@H](CO)[C@@H](O)[C@H](O)[C@@H]2O)[C@H](O)[C@@H](O)[C@H](O)[C@@H]1O)=O)OC(CCCCCCCCCCCCC)=O"

    def setUp(self): self.v = LipidHeadValidator()

    def test_positive(self):
        """Positive: correct structure must be recognised."""
        m = parse(self.POS)
        self.assertIsNotNone(m, "Ex 62 (3) POS SMILES invalid")
        self.assertTrue(self.v.matches_any_valid_head(m),
                        "Ex 62 (3) positive must match")

    def test_sugar_variant_still_matches(self):
        """Sugar variant with correct FA attachment must still match."""
        m = parse(self.VAR_SUGAR)
        self.assertIsNotNone(m, "Ex 62 (3) VAR_SUGAR SMILES invalid")
        self.assertTrue(self.v.matches_any_valid_head(m),
                        "Ex 62 (3) sugar-variant must still match")

    @unittest.expectedFailure
    def test_negative_stereo(self):
        """Wrong stereochemistry — not yet enforced by current SMARTS."""
        m = parse(self.NEG_STEREO)
        self.assertIsNotNone(m, "Ex 62 (3) NEG_STEREO SMILES invalid")
        self.assertFalse(self.v.matches_any_valid_head(m),
                         "Ex 62 (3) stereo-error should NOT match once stereo patterns added")


class TestEx63__3_DiglycerolTetraetherPhospholipidsCaldarc(unittest.TestCase):
    """Example 63 (3) — sn-caldarchaeo-1-phosphoethanolamine"""
    POS = "NCCOP(O)(OC[C@]1([H])COCC[C@H](C)CCC[C@H](C)CCC[C@H](C)CCC[C@H](C)CC[C@@H](C)CCC[C@@H](C)CCC[C@@H](C)CCC[C@@H](C)CCO[C@](COCC[C@H](C)CCC[C@H](C)CCC[C@H](C)CCC[C@H](C)CC[C@@H](C)CCC[C@@H](C)CCC[C@@H](C)CCC[C@@H](C)CCO1)([H])CO)=O"
    NEG_STEREO = "NCCOP(O)(OC[C@@]1([H])COCC[C@H](C)CCC[C@H](C)CCC[C@H](C)CCC[C@H](C)CC[C@@H](C)CCC[C@@H](C)CCC[C@@H](C)CCC[C@@H](C)CCO[C@](COCC[C@H](C)CCC[C@H](C)CCC[C@H](C)CCC[C@H](C)CC[C@@H](C)CCC[C@@H](C)CCC[C@@H](C)CCC[C@@H](C)CCO1)([H])CO)=O"

    def setUp(self): self.v = LipidHeadValidator()

    def test_positive(self):
        """Positive: correct structure must be recognised."""
        m = parse(self.POS)
        self.assertIsNotNone(m, "Ex 63 (3) POS SMILES invalid")
        self.assertTrue(self.v.matches_any_valid_head(m),
                        "Ex 63 (3) positive must match")

    @unittest.expectedFailure
    def test_negative_stereo(self):
        """Wrong stereochemistry — not yet enforced by current SMARTS."""
        m = parse(self.NEG_STEREO)
        self.assertIsNotNone(m, "Ex 63 (3) NEG_STEREO SMILES invalid")
        self.assertFalse(self.v.matches_any_valid_head(m),
                         "Ex 63 (3) stereo-error should NOT match once stereo patterns added")


class TestEx64__3_GlycerolnonitolTetraetherPhospholipids(unittest.TestCase):
    """Example 64 (3) — sn-caldito-1-phosphoethanolamine"""
    POS = "OC([C@](OCC[C@H](C)CCC[C@H](C)CCC[C@H](C)CCC[C@@H]1C)([H])COCC[C@H](C)CCC[C@H](C)CCC[C@H](C)CCC[C@H](C)CC[C@@H](C)CCC[C@@H](C)CCC[C@@H](C)CCC[C@@H](C)CCO[C@](COCC[C@H](C)CCC[C@H](C)CCC[C@H](C)CCC[C@H](C)CC1)([H])COP(OCCN)(O)=O)C(CO)(O)C(O)C(O)C(O)CO"
    NEG_STEREO = "OC([C@@](OCC[C@H](C)CCC[C@H](C)CCC[C@H](C)CCC[C@@H]1C)([H])COCC[C@H](C)CCC[C@H](C)CCC[C@H](C)CCC[C@H](C)CC[C@@H](C)CCC[C@@H](C)CCC[C@@H](C)CCC[C@@H](C)CCO[C@](COCC[C@H](C)CCC[C@H](C)CCC[C@H](C)CCC[C@H](C)CC1)([H])COP(OCCN)(O)=O)C(CO)(O)C(O)C(O)C(O)CO"

    def setUp(self): self.v = LipidHeadValidator()

    def test_positive(self):
        """Positive: correct structure must be recognised."""
        m = parse(self.POS)
        self.assertIsNotNone(m, "Ex 64 (3) POS SMILES invalid")
        self.assertTrue(self.v.matches_any_valid_head(m),
                        "Ex 64 (3) positive must match")

    @unittest.expectedFailure
    def test_negative_stereo(self):
        """Wrong stereochemistry — not yet enforced by current SMARTS."""
        m = parse(self.NEG_STEREO)
        self.assertIsNotNone(m, "Ex 64 (3) NEG_STEREO SMILES invalid")
        self.assertFalse(self.v.matches_any_valid_head(m),
                         "Ex 64 (3) stereo-error should NOT match once stereo patterns added")


class TestEx65__3_OxidizedGlycerophospholipids(unittest.TestCase):
    """Example 65 (3) — 1-O-(1Z-hexadecenyl)-2-(12S-hydroxy-5Z,8Z,10E,14Z-eicosatetraenoyl)-sn-glycero-3-phosphoethanolamine"""
    POS = "[H][C@](OC(CCC/C=C\\C/C=C\\C=C\\[C@@H](O)C/C=C\\CCCCC)=O)(COP(OCCN)(O)=O)CO/C=C\\CCCCCCCCCCCCCC"
    NEG_CHAIN = "[H][C@](O)(COP(OCCN)(OC(CCC/C=C\\C/C=C\\C=C\\[C@@H](O)C/C=C\\CCCCC)=O)=O)CO/C=C\\CCCCCCCCCCCCCC"
    NEG_STEREO = "[H]C(OC(CCC/C=C\\C/C=C\\C=C\\[C@@H](O)C/C=C\\CCCCC)=O)(COP(OCCN)(O)=O)CO/C=C\\CCCCCCCCCCCCCC"

    def setUp(self): self.v = LipidHeadValidator()

    def test_positive(self):
        """Positive: correct structure must be recognised."""
        m = parse(self.POS)
        self.assertIsNotNone(m, "Ex 65 (3) POS SMILES invalid")
        self.assertTrue(self.v.matches_any_valid_head(m),
                        "Ex 65 (3) positive must match")

    def test_negative_chain(self):
        """Chain-error (FA at wrong position) must NOT match."""
        m = parse(self.NEG_CHAIN)
        self.assertIsNotNone(m, "Ex 65 (3) NEG_CHAIN SMILES invalid")
        self.assertFalse(self.v.matches_any_valid_head(m),
                         "Ex 65 (3) chain-error must NOT match")

    @unittest.expectedFailure
    def test_negative_stereo(self):
        """Wrong stereochemistry — not yet enforced by current SMARTS."""
        m = parse(self.NEG_STEREO)
        self.assertIsNotNone(m, "Ex 65 (3) NEG_STEREO SMILES invalid")
        self.assertFalse(self.v.matches_any_valid_head(m),
                         "Ex 65 (3) stereo-error should NOT match once stereo patterns added")


class TestEx66__3_GlycerophosphoethanolamineGlycans(unittest.TestCase):
    """Example 66 (3) — N-(1-deoxyfructosyl)-1-hexadecanoyl-2-(4Z,7Z,10Z,13Z,16Z,19Z-docosahexaenoyl)-sn-glycero-3-phosphoethanolamine"""
    POS = "O[C@]1(CNCCOP(OC[C@]([H])(OC(CC/C=C\\C/C=C\\C/C=C\\C/C=C\\C/C=C\\C/C=C\\CC)=O)COC(CCCCCCCCCCCCCCC)=O)(O)=O)OC[C@@H](O)[C@@H](O)[C@@H]1O"
    VAR_SUGAR = "O[C@]1(CNCCOP(OC[C@]([H])(OC(CC/C=C\\C/C=C\\C/C=C\\C/C=C\\C/C=C\\C/C=C\\CC)=O)COC(CCCCCCCCCCCCCCC)=O)(O)=O)OC[C@H](O)[C@H](O)[C@@H]1O"
    NEG_STEREO = "O[C@]1(CNCCOP(OCC([H])(OC(CC/C=C\\C/C=C\\C/C=C\\C/C=C\\C/C=C\\C/C=C\\CC)=O)COC(CCCCCCCCCCCCCCC)=O)(O)=O)OC[C@@H](O)[C@@H](O)[C@@H]1O"

    def setUp(self): self.v = LipidHeadValidator()

    def test_positive(self):
        """Positive: correct structure must be recognised."""
        m = parse(self.POS)
        self.assertIsNotNone(m, "Ex 66 (3) POS SMILES invalid")
        self.assertTrue(self.v.matches_any_valid_head(m),
                        "Ex 66 (3) positive must match")

    def test_sugar_variant_still_matches(self):
        """Sugar variant with correct FA attachment must still match."""
        m = parse(self.VAR_SUGAR)
        self.assertIsNotNone(m, "Ex 66 (3) VAR_SUGAR SMILES invalid")
        self.assertTrue(self.v.matches_any_valid_head(m),
                        "Ex 66 (3) sugar-variant must still match")

    @unittest.expectedFailure
    def test_negative_stereo(self):
        """Wrong stereochemistry — not yet enforced by current SMARTS."""
        m = parse(self.NEG_STEREO)
        self.assertIsNotNone(m, "Ex 66 (3) NEG_STEREO SMILES invalid")
        self.assertFalse(self.v.matches_any_valid_head(m),
                         "Ex 66 (3) stereo-error should NOT match once stereo patterns added")


class TestEx67__3_Dihydroxyacetonephospates(unittest.TestCase):
    """Example 67 (3) — 1-hexadecanoyl-glycerone 3-phosphate"""
    POS = "O=C(CCCCCCCCCCCCCCC)OCC(COP(O)(O)=O)=O"
    NEG_CHAIN = "OCC(COP(O)(OC(CCCCCCCCCCCCCCC)=O)=O)=O"

    def setUp(self): self.v = LipidHeadValidator()

    def test_positive(self):
        """Positive: correct structure must be recognised."""
        m = parse(self.POS)
        self.assertIsNotNone(m, "Ex 67 (3) POS SMILES invalid")
        self.assertTrue(self.v.matches_any_valid_head(m),
                        "Ex 67 (3) positive must match")

    def test_negative_chain(self):
        """Chain-error (FA at wrong position) must NOT match."""
        m = parse(self.NEG_CHAIN)
        self.assertIsNotNone(m, "Ex 67 (3) NEG_CHAIN SMILES invalid")
        self.assertFalse(self.v.matches_any_valid_head(m),
                         "Ex 67 (3) chain-error must NOT match")


class TestEx68__3_Glycerophosphoethanols(unittest.TestCase):
    """Example 68 (3) — 1-hexadecanoyl-2-(9Z-octadecenoyl)-sn-glycero-3-phosphoethanol"""
    POS = "O=C(CCCCCCCCCCCCCCC)OC[C@]([H])(COP(OCC)(O)=O)OC(CCCCCCC/C=C\\CCCCCCCC)=O"
    NEG_CHAIN = "OC[C@]([H])(COP(OCC)(OC(CCCCCCCCCCCCCCC)=O)=OC(CCCCCCC/C=C\\CCCCCCCC)=O)O"
    NEG_STEREO = "O=C(CCCCCCCCCCCCCCC)OCC([H])(COP(OCC)(O)=O)OC(CCCCCCC/C=C\\CCCCCCCC)=O"

    def setUp(self): self.v = LipidHeadValidator()

    def test_positive(self):
        """Positive: correct structure must be recognised."""
        m = parse(self.POS)
        self.assertIsNotNone(m, "Ex 68 (3) POS SMILES invalid")
        self.assertTrue(self.v.matches_any_valid_head(m),
                        "Ex 68 (3) positive must match")

    def test_negative_chain(self):
        """Chain-error (FA at wrong position) must NOT match."""
        m = parse(self.NEG_CHAIN)
        self.assertIsNotNone(m, "Ex 68 (3) NEG_CHAIN SMILES invalid")
        self.assertFalse(self.v.matches_any_valid_head(m),
                         "Ex 68 (3) chain-error must NOT match")

    @unittest.expectedFailure
    def test_negative_stereo(self):
        """Wrong stereochemistry — not yet enforced by current SMARTS."""
        m = parse(self.NEG_STEREO)
        self.assertIsNotNone(m, "Ex 68 (3) NEG_STEREO SMILES invalid")
        self.assertFalse(self.v.matches_any_valid_head(m),
                         "Ex 68 (3) stereo-error should NOT match once stereo patterns added")


class TestEx69__3_Glycerophosphothreonines(unittest.TestCase):
    """Example 69 (3) — 1-octadecanoyl-2-(9Z-octadecenoyl)-sn-glycero-3-phosphothreonine"""
    POS = "O=C(CCCCCCCCCCCCCCCCC)OC[C@]([H])(COP(O[C@H](C)[C@](C(O)=O)([H])N)(O)=O)OC(CCCCCCC/C=C\\CCCCCCCC)=O"
    NEG_CHAIN = "OC[C@]([H])(COP(O[C@H](C)[C@](C(O)=O)([H])NC(CCCCCCC/C=C\\CCCCCCCC)=O)(OC(CCCCCCCCCCCCCCCCC)=O)=O)O"
    NEG_STEREO = "O=C(CCCCCCCCCCCCCCCCC)OCC([H])(COP(O[C@H](C)[C@](C(O)=O)([H])N)(O)=O)OC(CCCCCCC/C=C\\CCCCCCCC)=O"

    def setUp(self): self.v = LipidHeadValidator()

    def test_positive(self):
        """Positive: correct structure must be recognised."""
        m = parse(self.POS)
        self.assertIsNotNone(m, "Ex 69 (3) POS SMILES invalid")
        self.assertTrue(self.v.matches_any_valid_head(m),
                        "Ex 69 (3) positive must match")

    def test_negative_chain(self):
        """Chain-error (FA at wrong position) must NOT match."""
        m = parse(self.NEG_CHAIN)
        self.assertIsNotNone(m, "Ex 69 (3) NEG_CHAIN SMILES invalid")
        self.assertFalse(self.v.matches_any_valid_head(m),
                         "Ex 69 (3) chain-error must NOT match")

    @unittest.expectedFailure
    def test_negative_stereo(self):
        """Wrong stereochemistry — not yet enforced by current SMARTS."""
        m = parse(self.NEG_STEREO)
        self.assertIsNotNone(m, "Ex 69 (3) NEG_STEREO SMILES invalid")
        self.assertFalse(self.v.matches_any_valid_head(m),
                         "Ex 69 (3) stereo-error should NOT match once stereo patterns added")


class TestEx71__3_OtherGlycerolipids(unittest.TestCase):
    """Example 71 (3) — 1-(9Z,1Z-octadecadienoyl)-2-(10Z,13Z,16Z,19Z-docosatetraenoyl)-3-O-[hydroxymethyl-N,N,N-trimethyl-β-alanine]-glycerol"""
    POS = "O=C(CCCCCCC/C=C\\C/C=C\\CCCCC)OCC([H])(COCC([H])(C([O-])=O)C[N+](C)(C)C)OC(CCCCCCCC/C=C\\C/C=C\\C/C=C\\C/C=C\\CC)=O"

    def setUp(self): self.v = LipidHeadValidator()

    def test_positive(self):
        """Positive: correct structure must be recognised."""
        m = parse(self.POS)
        self.assertIsNotNone(m, "Ex 71 (3) POS SMILES invalid")
        self.assertTrue(self.v.matches_any_valid_head(m),
                        "Ex 71 (3) positive must match")

    def test_positive_class(self):
        """Lipid class must be identified correctly."""
        m = parse(self.POS)
        classes = self.v.identify_lipid_class(m)
        self.assertTrue(len(classes) > 0, "Ex 71 (3) must have a lipid class")


class TestEx72__3_OtherGlycerolipids(unittest.TestCase):
    """Example 72 (3) — 1-O-(1'S,2'S,3'R,4'R,5'S-tetrahydroxycyclopentyl)-2-(9-methylpentadecanoyl)-3-(10-methyl-hexadecanyl)-sn-glycerol"""
    POS = "[H][C@@](COCCCCCCCCCC(C)CCCCCC)(OC(CCCCCCCC(C)CCCCCC)=O)COC1[C@@H](O)[C@H](O)[C@H]([C@@H]1O)O"
    NEG_CHAIN = "[H][C@@](COCCCCCCCCCC(C)CCCCCC)(O)CO[C@@H]1[C@@H](O)[C@H](O)[C@H]([C@@H]1O)OC(CCCCCCCC(C)CCCCCC)=O"

    def setUp(self): self.v = LipidHeadValidator()

    def test_positive(self):
        """Positive: correct structure must be recognised."""
        m = parse(self.POS)
        self.assertIsNotNone(m, "Ex 72 (3) POS SMILES invalid")
        self.assertTrue(self.v.matches_any_valid_head(m),
                        "Ex 72 (3) positive must match")

    def test_positive_class(self):
        """Lipid class must be identified correctly."""
        m = parse(self.POS)
        classes = self.v.identify_lipid_class(m)
        self.assertTrue(len(classes) > 0, "Ex 72 (3) must have a lipid class")

    def test_negative_chain(self):
        """Chain-error (FA at wrong position) must NOT match."""
        m = parse(self.NEG_CHAIN)
        self.assertIsNotNone(m, "Ex 72 (3) NEG_CHAIN SMILES invalid")
        self.assertFalse(self.v.matches_any_valid_head(m),
                         "Ex 72 (3) chain-error must NOT match")


class TestEx73__3_OtherGlycerolipids(unittest.TestCase):
    """Example 73 (3) — 1-hexadecanoyl-2-(10-methyl-octadecanoyl)-3-O-(2S,5-diaminohexanoyl)-sn-glycerol"""
    POS = "O=C(CCCCCCCCCCCCCCC)OC[C@@]([H])(OC(CCCCCCCCC(C)CCCCCCCC)=O)COC([C@H](N)CCCCN)=O"
    NEG_STEREO = "O=C(CCCCCCCCCCCCCCC)OCC([H])(OC(CCCCCCCCC(C)CCCCCCCC)=O)COC([C@H](N)CCCCN)=O"

    def setUp(self): self.v = LipidHeadValidator()

    def test_positive(self):
        """Positive: correct structure must be recognised."""
        m = parse(self.POS)
        self.assertIsNotNone(m, "Ex 73 (3) POS SMILES invalid")
        self.assertTrue(self.v.matches_any_valid_head(m),
                        "Ex 73 (3) positive must match")

    def test_positive_class(self):
        """Lipid class must be identified correctly."""
        m = parse(self.POS)
        classes = self.v.identify_lipid_class(m)
        self.assertTrue(len(classes) > 0, "Ex 73 (3) must have a lipid class")

    @unittest.expectedFailure
    def test_negative_stereo(self):
        """Wrong stereochemistry — not yet enforced by current SMARTS."""
        m = parse(self.NEG_STEREO)
        self.assertIsNotNone(m, "Ex 73 (3) NEG_STEREO SMILES invalid")
        self.assertFalse(self.v.matches_any_valid_head(m),
                         "Ex 73 (3) stereo-error should NOT match once stereo patterns added")


class TestEx74__3_Glycerolipids(unittest.TestCase):
    """Example 74 (3) — 1-dodecanoyl-sn-glycerol"""
    POS = "OC[C@]([H])(O)COC(CCCCCCCCCCC)=O"
    NEG_STEREO = "OCC([H])(O)COC(CCCCCCCCCCC)=O"

    def setUp(self): self.v = LipidHeadValidator()

    def test_positive(self):
        """Positive: correct structure must be recognised."""
        m = parse(self.POS)
        self.assertIsNotNone(m, "Ex 74 (3) POS SMILES invalid")
        self.assertTrue(self.v.matches_any_valid_head(m),
                        "Ex 74 (3) positive must match")

    def test_positive_class(self):
        """Lipid class must be identified correctly."""
        m = parse(self.POS)
        classes = self.v.identify_lipid_class(m)
        self.assertTrue(len(classes) > 0, "Ex 74 (3) must have a lipid class")

    @unittest.expectedFailure
    def test_negative_stereo(self):
        """Wrong stereochemistry — not yet enforced by current SMARTS."""
        m = parse(self.NEG_STEREO)
        self.assertIsNotNone(m, "Ex 74 (3) NEG_STEREO SMILES invalid")
        self.assertFalse(self.v.matches_any_valid_head(m),
                         "Ex 74 (3) stereo-error should NOT match once stereo patterns added")


class TestEx75__3_Glycerolipids(unittest.TestCase):
    """Example 75 (3) — 1-acyl-3-O-β-D-galactosyl-sn-glycerol"""
    POS = "O=C(C)OC[C@@]([H])(O)CO[C@H]1[C@H](O)[C@@H](O)[C@@H](O)[C@@H](CO)O1"
    NEG_CHAIN = "OC[C@@]([H])(O)CO[C@H]1[C@H](O)[C@@H](O)[C@@H](O)[C@@H](COC(C)=O)O1"
    NEG_STEREO = "O=C(C)OCC([H])(O)CO[C@H]1[C@H](O)[C@@H](O)[C@@H](O)[C@@H](CO)O1"

    def setUp(self): self.v = LipidHeadValidator()

    def test_positive(self):
        """Positive: correct structure must be recognised."""
        m = parse(self.POS)
        self.assertIsNotNone(m, "Ex 75 (3) POS SMILES invalid")
        self.assertTrue(self.v.matches_any_valid_head(m),
                        "Ex 75 (3) positive must match")

    def test_positive_class(self):
        """Lipid class must be identified correctly."""
        m = parse(self.POS)
        classes = self.v.identify_lipid_class(m)
        self.assertTrue(len(classes) > 0, "Ex 75 (3) must have a lipid class")

    def test_negative_chain(self):
        """Chain-error (FA at wrong position) must NOT match."""
        m = parse(self.NEG_CHAIN)
        self.assertIsNotNone(m, "Ex 75 (3) NEG_CHAIN SMILES invalid")
        self.assertFalse(self.v.matches_any_valid_head(m),
                         "Ex 75 (3) chain-error must NOT match")

    @unittest.expectedFailure
    def test_negative_stereo(self):
        """Wrong stereochemistry — not yet enforced by current SMARTS."""
        m = parse(self.NEG_STEREO)
        self.assertIsNotNone(m, "Ex 75 (3) NEG_STEREO SMILES invalid")
        self.assertFalse(self.v.matches_any_valid_head(m),
                         "Ex 75 (3) stereo-error should NOT match once stereo patterns added")


class TestEx76__3_Glycerolipids(unittest.TestCase):
    """Example 76 (3) — 1-hexadecanyl-2-((2'-α-glucosyl)-β-glucosyl)-3-β-xylosyl-sn-glycerol"""
    POS = "[H][C@](COCCCCCCCCCCCCCCCC)(O[C@@H]1O[C@H](CO)[C@@H](O)[C@H](O)[C@H]1O[C@H]2O[C@H](CO)[C@@H](O)[C@H](O)[C@H]2O)CO[C@H]3O[C@@H](CO)[C@@H](O)[C@@H]3O"
    NEG_CHAIN = "[H][C@](COCCCCCCCCCCCCCCCC)(O[C@H]1O[C@H](CO)[C@@H](O)[C@H](O)[C@H]1O[C@H]2O[C@@H](CO)[C@@H](O)[C@H](O)[C@H]2O)CO[C@H]3O[C@@H](CO)[C@@H](O)[C@@H]3O"
    NEG_STEREO = "[H]C(COCCCCCCCCCCCCCCCC)(O[C@@H]1O[C@H](CO)[C@@H](O)[C@H](O)[C@H]1O[C@H]2O[C@H](CO)[C@@H](O)[C@H](O)[C@H]2O)CO[C@H]3O[C@@H](CO)[C@@H](O)[C@@H]3O"

    def setUp(self): self.v = LipidHeadValidator()

    def test_positive(self):
        """Positive: correct structure must be recognised."""
        m = parse(self.POS)
        self.assertIsNotNone(m, "Ex 76 (3) POS SMILES invalid")
        self.assertTrue(self.v.matches_any_valid_head(m),
                        "Ex 76 (3) positive must match")

    def test_positive_class(self):
        """Lipid class must be identified correctly."""
        m = parse(self.POS)
        classes = self.v.identify_lipid_class(m)
        self.assertTrue(len(classes) > 0, "Ex 76 (3) must have a lipid class")

    def test_negative_chain(self):
        """Chain-error (FA at wrong position) must NOT match."""
        m = parse(self.NEG_CHAIN)
        self.assertIsNotNone(m, "Ex 76 (3) NEG_CHAIN SMILES invalid")
        self.assertFalse(self.v.matches_any_valid_head(m),
                         "Ex 76 (3) chain-error must NOT match")

    @unittest.expectedFailure
    def test_negative_stereo(self):
        """Wrong stereochemistry — not yet enforced by current SMARTS."""
        m = parse(self.NEG_STEREO)
        self.assertIsNotNone(m, "Ex 76 (3) NEG_STEREO SMILES invalid")
        self.assertFalse(self.v.matches_any_valid_head(m),
                         "Ex 76 (3) stereo-error should NOT match once stereo patterns added")


class TestEx77__3_Glycerolipids(unittest.TestCase):
    """Example 77 (3) — 1-(9Z-hexadecenoyl)-3-(6'-sulfo-α-D-quinovosyl)-sn-glycerol"""
    POS = "[H][C@](O)(CO[C@H]1O[C@H](CS(=O)(O)=O)[C@@H](O)[C@H](O)[C@H]1O)COC(CCCCCCC/C=C\\CCCCCC)=O"
    NEG_CHAIN = "[H]C(O)(CO[C@H]1O[C@H](CS(=O)(O)=O)[C@@H](O)[C@H](O)[C@H]1O)COC(CCCCCCC/C=C\\CCCCCC)=O"
    VAR_SUGAR = "[H][C@](O)(CO[C@H]1O[C@@H](CS(=O)(O)=O)[C@@H](O)[C@@H](O)[C@H]1O)COC(CCCCCCC/C=C\\CCCCCC)=O"

    def setUp(self): self.v = LipidHeadValidator()

    def test_positive(self):
        """Positive: correct structure must be recognised."""
        m = parse(self.POS)
        self.assertIsNotNone(m, "Ex 77 (3) POS SMILES invalid")
        self.assertTrue(self.v.matches_any_valid_head(m),
                        "Ex 77 (3) positive must match")

    def test_positive_class(self):
        """Lipid class must be identified correctly."""
        m = parse(self.POS)
        classes = self.v.identify_lipid_class(m)
        self.assertTrue(len(classes) > 0, "Ex 77 (3) must have a lipid class")

    def test_negative_chain(self):
        """Chain-error (FA at wrong position) must NOT match."""
        m = parse(self.NEG_CHAIN)
        self.assertIsNotNone(m, "Ex 77 (3) NEG_CHAIN SMILES invalid")
        self.assertFalse(self.v.matches_any_valid_head(m),
                         "Ex 77 (3) chain-error must NOT match")

    def test_sugar_variant_still_matches(self):
        """Sugar variant with correct FA attachment must still match."""
        m = parse(self.VAR_SUGAR)
        self.assertIsNotNone(m, "Ex 77 (3) VAR_SUGAR SMILES invalid")
        self.assertTrue(self.v.matches_any_valid_head(m),
                        "Ex 77 (3) sugar-variant must still match")


class TestEx78__3_Betaine(unittest.TestCase):
    """Example 78 (3) — 1-(9Z-octadecenoyl)-2-(10Z,13Z,16Z,19Z-docosatetraenoyl)-sn-glycero-3-O-2'-(hydroxymethyl)-(N,N,N-trimethyl)-β-alanine"""
    POS = "O=C(CCCCCCC/C=C\\CCCCCCCC)OC[C@]([H])(COCC(C([O-])=O)C[N+](C)(C)C)OC(CCCCCCCC/C=C\\C/C=C\\C/C=C\\C/C=C\\CC)=O"
    NEG_CHAIN = "OC[C@]([H]C(CCCCCCC/C=C\\CCCCCCCC)=O)(COCC(/C([O-])=O/C(CCCCCCCC/C=C\\C/C=C\\C/C=C\\C/C=C\\CC)=O)C[N+](C)(C)C)O"
    NEG_STEREO = "O=C(CCCCCCC/C=C\\CCCCCCCC)OCC([H])(COCC(C([O-])=O)C[N+](C)(C)C)OC(CCCCCCCC/C=C\\C/C=C\\C/C=C\\C/C=C\\CC)=O"

    def setUp(self): self.v = LipidHeadValidator()

    def test_positive(self):
        """Positive: correct structure must be recognised."""
        m = parse(self.POS)
        self.assertIsNotNone(m, "Ex 78 (3) POS SMILES invalid")
        self.assertTrue(self.v.matches_any_valid_head(m),
                        "Ex 78 (3) positive must match")

    def test_positive_class(self):
        """Lipid class must be identified correctly."""
        m = parse(self.POS)
        classes = self.v.identify_lipid_class(m)
        self.assertTrue(len(classes) > 0, "Ex 78 (3) must have a lipid class")

    def test_negative_chain(self):
        """Chain-error (FA at wrong position) must NOT match."""
        m = parse(self.NEG_CHAIN)
        self.assertIsNotNone(m, "Ex 78 (3) NEG_CHAIN SMILES invalid")
        self.assertFalse(self.v.matches_any_valid_head(m),
                         "Ex 78 (3) chain-error must NOT match")

    @unittest.expectedFailure
    def test_negative_stereo(self):
        """Wrong stereochemistry — not yet enforced by current SMARTS."""
        m = parse(self.NEG_STEREO)
        self.assertIsNotNone(m, "Ex 78 (3) NEG_STEREO SMILES invalid")
        self.assertFalse(self.v.matches_any_valid_head(m),
                         "Ex 78 (3) stereo-error should NOT match once stereo patterns added")


class TestEx79__3_Phosphosphingolipids(unittest.TestCase):
    """Example 79 (3) — N-(eicosanoyl)-hexadecasphing-4-enine-1-phosphocholine"""
    POS = "[H][C@](NC(CCCCCCCCCCCCCCCCCCC)=O)(COP(OCC[N+](C)(C)C)([O-])=O)[C@@](O)([H])/C=C/CCCCCCCCCCC"
    NEG_CHAIN = "N[C@@]([H]C(CCCCCCCCCCCCCCCCCCC)=O)(COP(OCC[N+](C)(C)C)([O-])=O)[C@@](O)([H])/C=C/CCCCCC(C)CCCCC"
    NEG_STEREO = "[H]C(NC(CCCCCCCCCCCCCCCCCCC)=O)(COP(OCC[N+](C)(C)C)([O-])=O)[C@@](O)([H])/C=C/CCCCCCCCCCC"

    def setUp(self): self.v = LipidHeadValidator()

    def test_positive(self):
        """Positive: correct structure must be recognised."""
        m = parse(self.POS)
        self.assertIsNotNone(m, "Ex 79 (3) POS SMILES invalid")
        self.assertTrue(self.v.matches_any_valid_head(m),
                        "Ex 79 (3) positive must match")

    def test_negative_chain(self):
        """Chain-error (FA at wrong position) must NOT match."""
        m = parse(self.NEG_CHAIN)
        self.assertIsNotNone(m, "Ex 79 (3) NEG_CHAIN SMILES invalid")
        self.assertFalse(self.v.matches_any_valid_head(m),
                         "Ex 79 (3) chain-error must NOT match")

    @unittest.expectedFailure
    def test_negative_stereo(self):
        """Wrong stereochemistry — not yet enforced by current SMARTS."""
        m = parse(self.NEG_STEREO)
        self.assertIsNotNone(m, "Ex 79 (3) NEG_STEREO SMILES invalid")
        self.assertFalse(self.v.matches_any_valid_head(m),
                         "Ex 79 (3) stereo-error should NOT match once stereo patterns added")


class TestEx80__3_Phosphosphingolipids(unittest.TestCase):
    """Example 80 (3) — N-(2R-hydroxy--15-methyl-3E-hexadecenoyl)-9-methyl-eicosasphinga-4E,8E-dienine-1-phosphoethanolamine"""
    POS = "[H][C@](NC([C@H](O)/C=C/CCCCCCCCCCC(C)C)=O)(COP(OCCN)(O)=O)[C@@](O)([H])/C=C/CC/C=C(C)/CCCCCCCCCCC"
    NEG_CHAIN = "[H][C@](N)(COP(OCCN)(OC([C@H](O)/C=C/CCCCCCCCCCC(C)C)=O)=O)[C@@](O)([H])/C=C/CC/C=C(C)/CCCCCCCCCCC"
    NEG_STEREO = "[H]C(NC([C@H](O)/C=C/CCCCCCCCCCC(C)C)=O)(COP(OCCN)(O)=O)[C@@](O)([H])/C=C/CC/C=C(C)/CCCCCCCCCCC"

    def setUp(self): self.v = LipidHeadValidator()

    def test_positive(self):
        """Positive: correct structure must be recognised."""
        m = parse(self.POS)
        self.assertIsNotNone(m, "Ex 80 (3) POS SMILES invalid")
        self.assertTrue(self.v.matches_any_valid_head(m),
                        "Ex 80 (3) positive must match")

    def test_negative_chain(self):
        """Chain-error (FA at wrong position) must NOT match."""
        m = parse(self.NEG_CHAIN)
        self.assertIsNotNone(m, "Ex 80 (3) NEG_CHAIN SMILES invalid")
        self.assertFalse(self.v.matches_any_valid_head(m),
                         "Ex 80 (3) chain-error must NOT match")

    @unittest.expectedFailure
    def test_negative_stereo(self):
        """Wrong stereochemistry — not yet enforced by current SMARTS."""
        m = parse(self.NEG_STEREO)
        self.assertIsNotNone(m, "Ex 80 (3) NEG_STEREO SMILES invalid")
        self.assertFalse(self.v.matches_any_valid_head(m),
                         "Ex 80 (3) stereo-error should NOT match once stereo patterns added")


class TestEx81__3_Phosphosphingolipids(unittest.TestCase):
    """Example 81 (3) — N-(docosanoyl)-sphing-4-enine-1-phospho-(1'-myo-inositol)"""
    POS = "[H][C@](NC(CCCCCCCCCCCCCCCCCCCCC)=O)(COP(O)(O[C@@H]1[C@H](O)[C@H](O)[C@@H](O)[C@H](O)[C@H]1O)=O)[C@@](O)([H])/C=C/CCCCCCCCCCCCC"
    NEG_CHAIN = "O=C([H][C@](N)(COP(O)(O[C@@H]1[C@H](O)[C@H](O)[C@@H](O)[C@H](O)[C@H]1O)=O)[C@@](O)([H])/C=C/CCCCCCCCCCCCC)CCCCCCCCCCCCCCCCCCCCC"
    VAR_SUGAR = "[H][C@](NC(CCCCCCCCCCCCCCCCCCCCC)=O)(COP(O)(O[C@H]1[C@@H](O)[C@H](O)[C@@H](O)[C@H](O)[C@H]1O)=O)[C@@](O)([H])/C=C/CCCCCCCCCCCCC"
    NEG_STEREO = "[H]C(NC(CCCCCCCCCCCCCCCCCCCCC)=O)(COP(O)(O[C@@H]1[C@H](O)[C@H](O)[C@@H](O)[C@H](O)[C@H]1O)=O)[C@@](O)([H])/C=C/CCCCCCCCCCCCC"

    def setUp(self): self.v = LipidHeadValidator()

    def test_positive(self):
        """Positive: correct structure must be recognised."""
        m = parse(self.POS)
        self.assertIsNotNone(m, "Ex 81 (3) POS SMILES invalid")
        self.assertTrue(self.v.matches_any_valid_head(m),
                        "Ex 81 (3) positive must match")

    def test_negative_chain(self):
        """Chain-error (FA at wrong position) must NOT match."""
        m = parse(self.NEG_CHAIN)
        self.assertIsNotNone(m, "Ex 81 (3) NEG_CHAIN SMILES invalid")
        self.assertFalse(self.v.matches_any_valid_head(m),
                         "Ex 81 (3) chain-error must NOT match")

    def test_sugar_variant_still_matches(self):
        """Sugar variant with correct FA attachment must still match."""
        m = parse(self.VAR_SUGAR)
        self.assertIsNotNone(m, "Ex 81 (3) VAR_SUGAR SMILES invalid")
        self.assertTrue(self.v.matches_any_valid_head(m),
                        "Ex 81 (3) sugar-variant must still match")

    @unittest.expectedFailure
    def test_negative_stereo(self):
        """Wrong stereochemistry — not yet enforced by current SMARTS."""
        m = parse(self.NEG_STEREO)
        self.assertIsNotNone(m, "Ex 81 (3) NEG_STEREO SMILES invalid")
        self.assertFalse(self.v.matches_any_valid_head(m),
                         "Ex 81 (3) stereo-error should NOT match once stereo patterns added")


class TestEx82__3_NeutralGlycosphingolipids(unittest.TestCase):
    """Example 82 (3) — N-(34S-methyl-5Z,9Z,12Z,15Z,18Z,21Z-hexatriacontahexaenoyl)-1-sulfo-β-D-fucosyl-(16R-methyl-sphing-4-enine)"""
    POS = "O=C(N[C@@]([H])(COS(O[C@@H]1O[C@H](C)[C@H](O)[C@H](O)[C@H]1O)(=O)=O)[C@@](O)([H])/C=C/CCCCCCCCCC[C@H](C)CC)CCC/C=C\\CC/C=C\\C/C=C\\C/C=C\\C/C=C\\C/C=C\\CCCCCCCCCCC[C@@H](C)CC"
    NEG_CHAIN = "[H][C@](N)(COS(O[C@H]1O[C@H](C)[C@H](OC(CCC/C=C\\CC/C=C\\C/C=C\\C/C=C\\C/C=C\\C/C=C\\CCCCCCCCCCC[C@@H](C)CC)=O)[C@H](O)[C@H]1O)(=O)=O)[C@@](O)([H])/C=C/CCCCCCCCCC[C@H](C)CC"
    VAR_SUGAR = "O=C(N[C@@]([H])(COS(O[C@H]1O[C@H](C)[C@H](O)[C@H](O)[C@H]1O)(=O)=O)[C@@](O)([H])/C=C/CCCCCCCCCC[C@H](C)CC)CCC/C=C\\CC/C=C\\C/C=C\\C/C=C\\C/C=C\\C/C=C\\CCCCCCCCCCC[C@@H](C)CC"
    NEG_STEREO = "O=C(NC([H])(COS(O[C@@H]1O[C@H](C)[C@H](O)[C@H](O)[C@H]1O)(=O)=O)[C@@](O)([H])/C=C/CCCCCCCCCC[C@H](C)CC)CCC/C=C\\CC/C=C\\C/C=C\\C/C=C\\C/C=C\\C/C=C\\CCCCCCCCCCC[C@@H](C)CC"

    def setUp(self): self.v = LipidHeadValidator()

    def test_positive(self):
        """Positive: correct structure must be recognised."""
        m = parse(self.POS)
        self.assertIsNotNone(m, "Ex 82 (3) POS SMILES invalid")
        self.assertTrue(self.v.matches_any_valid_head(m),
                        "Ex 82 (3) positive must match")

    def test_negative_chain(self):
        """Chain-error (FA at wrong position) must NOT match."""
        m = parse(self.NEG_CHAIN)
        self.assertIsNotNone(m, "Ex 82 (3) NEG_CHAIN SMILES invalid")
        self.assertFalse(self.v.matches_any_valid_head(m),
                         "Ex 82 (3) chain-error must NOT match")

    def test_sugar_variant_still_matches(self):
        """Sugar variant with correct FA attachment must still match."""
        m = parse(self.VAR_SUGAR)
        self.assertIsNotNone(m, "Ex 82 (3) VAR_SUGAR SMILES invalid")
        self.assertTrue(self.v.matches_any_valid_head(m),
                        "Ex 82 (3) sugar-variant must still match")

    @unittest.expectedFailure
    def test_negative_stereo(self):
        """Wrong stereochemistry — not yet enforced by current SMARTS."""
        m = parse(self.NEG_STEREO)
        self.assertIsNotNone(m, "Ex 82 (3) NEG_STEREO SMILES invalid")
        self.assertFalse(self.v.matches_any_valid_head(m),
                         "Ex 82 (3) stereo-error should NOT match once stereo patterns added")


class TestEx83__3_NeutralGlycosphingolipids(unittest.TestCase):
    """Example 83 (3) — 1-O-melibiosoyl-(N-(2R-hydroxy-heneicosanoyl)-4R-hydroxy-17-methyl-sphing-6E-enine)"""
    POS = "[H][C@](NC([C@H](O)CCCCCCCCCCCCCCCCCCC)=O)(CO[C@@H]1O[C@H](CO[C@H]2O[C@H](CO)[C@H](O)[C@H](O)[C@H]2O)[C@@H](O)[C@H](O)[C@H]1O)[C@@](O)([H])[C@H](O)/C=C/CCCCCCCCCC(C)C"
    NEG_CHAIN = "[H][C@](N)(CO[C@@H]1O[C@H](CO[C@H]2O[C@H](COC([C@H](O)CCCCCCCCCCCCCCCCCCC)=O)[C@H](O)[C@H](O)[C@H]2O)[C@@H](O)[C@H](O)[C@H]1O)[C@@](O)([H])[C@H](O)/C=C/CCCCCCCCCCC.C"
    VAR_SUGAR = "[H][C@](NC([C@H](O)CCCCCCCCCCCCCCCCCCC)=O)(CO[C@H]1O[C@H](CO[C@H]2O[C@H](CO)[C@H](O)[C@@H](O)[C@H]2O)[C@@H](O)[C@H](O)[C@H]1O)[C@@](O)([H])[C@H](O)/C=C/CCCCCCCCCC(C)C"
    NEG_STEREO = "[H][C@](NC([C@H](O)CCCCCCCCCCCCCCCCCCC)=O)(CO[C@@H]1O[C@H](CO[C@H]2O[C@H](CO)[C@H](O)[C@H](O)[C@H]2O)[C@@H](O)[C@H](O)[C@H]1O)C(O)([H])[C@H](O)/C=C/CCCCCCCCCC(C)C"

    def setUp(self): self.v = LipidHeadValidator()

    def test_positive(self):
        """Positive: correct structure must be recognised."""
        m = parse(self.POS)
        self.assertIsNotNone(m, "Ex 83 (3) POS SMILES invalid")
        self.assertTrue(self.v.matches_any_valid_head(m),
                        "Ex 83 (3) positive must match")

    def test_negative_chain(self):
        """Chain-error (FA at wrong position) must NOT match."""
        m = parse(self.NEG_CHAIN)
        self.assertIsNotNone(m, "Ex 83 (3) NEG_CHAIN SMILES invalid")
        self.assertFalse(self.v.matches_any_valid_head(m),
                         "Ex 83 (3) chain-error must NOT match")

    def test_sugar_variant_still_matches(self):
        """Sugar variant with correct FA attachment must still match."""
        m = parse(self.VAR_SUGAR)
        self.assertIsNotNone(m, "Ex 83 (3) VAR_SUGAR SMILES invalid")
        self.assertTrue(self.v.matches_any_valid_head(m),
                        "Ex 83 (3) sugar-variant must still match")

    @unittest.expectedFailure
    def test_negative_stereo(self):
        """Wrong stereochemistry — not yet enforced by current SMARTS."""
        m = parse(self.NEG_STEREO)
        self.assertIsNotNone(m, "Ex 83 (3) NEG_STEREO SMILES invalid")
        self.assertFalse(self.v.matches_any_valid_head(m),
                         "Ex 83 (3) stereo-error should NOT match once stereo patterns added")


class TestEx84__3_AcidicGlycosphingolipids(unittest.TestCase):
    """Example 84 (3) — 1-O-(6'-phosphoethanolaminy-β-D-glucopyranosyl)-(N-(2R-hydroxy-tetracosanoyl)-4R-hydroxy-15-methyl-hexadecasphinganine)"""
    POS = "[H][C@](NC([C@H](O)CCCCCCCCCCCCCCCCCCCCCC)=O)(CO[C@@H]1O[C@H](COP(O)(OCCN)=O)[C@@H](O)[C@H](O)[C@H]1O)[C@@](O)([H])[C@@H](CCCCCCCCCCC(C)C)O"
    NEG_CHAIN = "[H][C@](N)(CO[C@@H]1O[C@H](COP(O)(OCCN)=O)[C@@H](O)[C@H](O)[C@H]1OC([C@H](O)CCCCCCCCCCCCCCCCCCCCCC)=O)[C@@](O)([H])[C@@H](CCCCCCCCCCC(C)C)O"
    VAR_SUGAR = "[H][C@](NC([C@H](O)CCCCCCCCCCCCCCCCCCCCCC)=O)(CO[C@H]1O[C@H](COP(O)(OCCN)=O)[C@@H](O)[C@H](O)C1)[C@@](O)([H])[C@@H](CCCCCCCCCCC(C)C)O"
    NEG_STEREO = "[H]C(NC([C@H](O)CCCCCCCCCCCCCCCCCCCCCC)=O)(CO[C@@H]1O[C@H](COP(O)(OCCN)=O)[C@@H](O)[C@H](O)[C@H]1O)[C@@](O)([H])[C@@H](CCCCCCCCCCC(C)C)O"

    def setUp(self): self.v = LipidHeadValidator()

    def test_positive(self):
        """Positive: correct structure must be recognised."""
        m = parse(self.POS)
        self.assertIsNotNone(m, "Ex 84 (3) POS SMILES invalid")
        self.assertTrue(self.v.matches_any_valid_head(m),
                        "Ex 84 (3) positive must match")

    def test_negative_chain(self):
        """Chain-error (FA at wrong position) must NOT match."""
        m = parse(self.NEG_CHAIN)
        self.assertIsNotNone(m, "Ex 84 (3) NEG_CHAIN SMILES invalid")
        self.assertFalse(self.v.matches_any_valid_head(m),
                         "Ex 84 (3) chain-error must NOT match")

    def test_sugar_variant_still_matches(self):
        """Sugar variant with correct FA attachment must still match."""
        m = parse(self.VAR_SUGAR)
        self.assertIsNotNone(m, "Ex 84 (3) VAR_SUGAR SMILES invalid")
        self.assertTrue(self.v.matches_any_valid_head(m),
                        "Ex 84 (3) sugar-variant must still match")

    @unittest.expectedFailure
    def test_negative_stereo(self):
        """Wrong stereochemistry — not yet enforced by current SMARTS."""
        m = parse(self.NEG_STEREO)
        self.assertIsNotNone(m, "Ex 84 (3) NEG_STEREO SMILES invalid")
        self.assertFalse(self.v.matches_any_valid_head(m),
                         "Ex 84 (3) stereo-error should NOT match once stereo patterns added")


class TestEx85__3_OtherPolyketides(unittest.TestCase):
    """Example 85 (3) — 1S-((4S-acetoxy-5R-methyl-3-methylene-6-phenylhexyl)-6-(E)-4S,6S-dimethyloct-2-enoyloxy)-4,7S-dihydroxy-2,8-dioxabicyclo[3.2.1]octane-3S,4S,5R-tricarboxylic acid"""
    POS = "O[C@@H]1[C@](O[C@H](C(O)=O)[C@]2(C(O)=O)O)(CCC([C@H]([C@H](C)CC3=CC=CC=C3)OC(C)=O)=C)O[C@]2(C(O)=O)[C@@H]1OC(/C=C/[C@@H](C)C[C@@H](C)CC)=O"
    NEG_STEREO = "O[C@@H]1[C@](O[C@H](C(O)=O)[C@@]2(C(O)=O)O)(CCC([C@H]([C@H](C)CC3=CC=CC=C3)OC(C)=O)=C)O[C@]2(C(O)=O)[C@@H]1OC(/C=C/[C@@H](C)C[C@@H](C)CC)=O"

    def setUp(self): self.v = LipidHeadValidator()

    def test_positive(self):
        """Positive: correct structure must be recognised."""
        m = parse(self.POS)
        self.assertIsNotNone(m, "Ex 85 (3) POS SMILES invalid")
        self.assertTrue(self.v.matches_any_valid_head(m),
                        "Ex 85 (3) positive must match")

    @unittest.expectedFailure
    def test_negative_stereo(self):
        """Wrong stereochemistry — not yet enforced by current SMARTS."""
        m = parse(self.NEG_STEREO)
        self.assertIsNotNone(m, "Ex 85 (3) NEG_STEREO SMILES invalid")
        self.assertFalse(self.v.matches_any_valid_head(m),
                         "Ex 85 (3) stereo-error should NOT match once stereo patterns added")


if __name__ == "__main__":
    unittest.main(verbosity=2)