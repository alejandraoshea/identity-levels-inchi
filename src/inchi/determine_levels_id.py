from rdkit import Chem
from rdkit.Chem import rdFMCS
from rdkit.Chem.SaltRemover import SaltRemover
from inchi.inchi_parser import InChiParser
from inchi.inchi_layers_enum import InchiLayers

class InChi:
    def isCompleteIdentity(inchi1: str, inchi2: str) -> bool:
        return (
            InChi.main_layer(inchi1, inchi2) and
            InChi.charge_layer(inchi1, inchi2) and
            InChi.stereochemical_layer(inchi1, inchi2) and
            InChi.isotopic_layer(inchi1, inchi2) and
            InChi.fixed_H_layer(inchi1, inchi2) and
            InChi.reconnected_layer(inchi1, inchi2)
        )

    #main layer (sublayers: atom connections and hydrogen atoms)
    def main_layer(inchi1: str, inchi2: str) -> bool:
        return(
            InChiParser.getMainLayer(inchi1) == InChiParser.getMainLayer(inchi2) and
            InChiParser.getAtomConnectionsSublayer(inchi1) == InChiParser.getAtomConnectionsSublayer(inchi2) and
            InChiParser.getHydrogenAtomsSublayer(inchi1) == InChiParser.getHydrogenAtomsSublayer(inchi2)
        )
    
    #charge layer (sublayers: charge and proton)
    def charge_layer(inchi1: str, inchi2: str) -> bool:
        return (
            InChiParser.getChargeSublayer(inchi1) == InChiParser.getChargeSublayer(inchi2) and
            InChiParser.getProtonSublayer(inchi1) == InChiParser.getProtonSublayer(inchi2)
        )
    
    def compare_charge_effects(original_inchi, neutral_inchi):
        return {
            "orig_no_charge": InChiParser.removeChargeLayers(original_inchi),
            "neutral_no_charge": InChiParser.removeChargeLayers(neutral_inchi),
            "charge_independent_equal":
                InChiParser.removeChargeLayers(original_inchi)
                == InChiParser.removeChargeLayers(neutral_inchi)
        }
    
    #stereochemical layer (sublayers: double bonds, tetrahedrals, type)
    def stereochemical_layer(inchi1: str, inchi2: str) -> bool:
        return (
            InChiParser.getDoubleBondsSublayer(inchi1) == InChiParser.getDoubleBondsSublayer(inchi2) and
            InChiParser.getTetrahedralStereoSublayer(inchi1) == InChiParser.getTetrahedralStereoSublayer(inchi2) and
            InChiParser.getTypeStereoInfoSublayer(inchi1) == InChiParser.getTypeStereoInfoSublayer(inchi2)
        )
    
    def isotopic_layer(inchi1: str, inchi2: str) -> bool:
        return (
            InChiParser.getIsotopicLayer(inchi1) == InChiParser.getIsotopicLayer(inchi2) and
            InChiParser.getIsotopicHydrogenSublayer(inchi1) == InChiParser.getIsotopicHydrogenSublayer(inchi2) and
            InChiParser.getIsotopicStereoSublayer(inchi1) == InChiParser.getIsotopicStereoSublayer(inchi2)
        )
    #Isotopic layer (prefix: "i"), may include sublayers:[13]

    def fixed_H_layer(inchi1: str, inchi2: str) -> bool:
        return InChiParser.getFixedHLayer(inchi1) == InChiParser.getFixedHLayer(inchi2)

    def reconnected_layer(inchi1: str, inchi2: str) -> bool:
        return InChiParser.getReconnectedLayer(inchi1) == InChiParser.getReconnectedLayer(inchi2)
    # never included in standard InChI


    def mol_from_inchi(inchi: str):
        try:
            mol = Chem.MolFromInchi(inchi)
            if mol is None:
                raise ValueError(f"Invalid InChI: {inchi}")
            return mol
        except Exception as e:
            print(f"RDKit conversion failed for InChI: {inchi}\nError: {e}")
            return None
    
    def main_fragment(mol):
        remover = SaltRemover()
        try:
            mol_clean = remover.StripMol(mol, dontRemoveEverything=True)
            Chem.AssignStereochemistry(mol_clean, cleanIt=True, force=True)
            return mol_clean
        except Exception:
            return mol
    
    def areEqualDisolvedSalts(inchi1: str, inchi2: str) -> bool:
        mol1 = InChi.mol_from_inchi(inchi1)
        mol2 = InChi.mol_from_inchi(inchi2)
        if not mol1 or not mol2:
            return False    
        main1 = InChi.main_fragment(mol1)
        main2 = InChi.main_fragment(mol2)      
        return Chem.MolToInchi(main1) == Chem.MolToInchi(main2)

    """def areEqualNoCharges(inchi1: str, inchi2:str) -> bool:
        inchi1_no_charge = InChiParser.removeChargeLayers(inchi1)
        inchi2_no_charge = InChiParser.removeChargeLayers(inchi2)
        return inchi1_no_charge == inchi2_no_charge
        """
    
    @staticmethod
    def compare_after_neutralization(inchi1: str, inchi2: str) -> bool:
        mol1 = InChi.mol_from_inchi(inchi1)
        mol2 = InChi.mol_from_inchi(inchi2)

        if mol1 is None or mol2 is None:
            return False

        neutral1 = InChiParser.neutralize_molecule(mol1)
        neutral2 = InChiParser.neutralize_molecule(mol2)

        # we compare connectivity-only signatures
        sig1 = Chem.MolToSmiles(neutral1, canonical=True, isomericSmiles=False,)
        sig2 = Chem.MolToSmiles(neutral2,canonical=True,isomericSmiles=False,)
        
        return sig1 == sig2

    def areEqualNoCharges(inchi1: str, inchi2: str) -> bool:
        if inchi1 is None or inchi2 is None:
            return False

        # remove salts
        mol1 = InChi.mol_from_inchi(inchi1)
        mol2 = InChi.mol_from_inchi(inchi2)

        if mol1 and mol2:
            mol1 = InChi.main_fragment(mol1)
            mol2 = InChi.main_fragment(mol2)
            inchi1 = Chem.MolToInchi(mol1)
            inchi2 = Chem.MolToInchi(mol2)

        # remove isotopes
        inchi1 = InChiParser.removeIsotopicLayers(inchi1)
        inchi2 = InChiParser.removeIsotopicLayers(inchi2)

        # remove charge layer
        inchi1 = InChiParser.removeChargeLayers(inchi1)
        inchi2 = InChiParser.removeChargeLayers(inchi2)

        return inchi1 == inchi2

    def areEqualNoStereo(inchi1: str, inchi2: str) -> bool:
        inchi1_no_stereo = InChiParser.removeStereoLayers(inchi1)
        inchi2_no_stereo = InChiParser.removeStereoLayers(inchi2)
        return inchi1_no_stereo == inchi2_no_stereo

    #stereochemical layer - sublayer
    def areEqualNoPositionDoubleBond(inchi1: str, inchi2:str) -> bool:
        inchi1_no_double_bonds = InChiParser.removeDoubleBondsSublayer(inchi1)
        inchi2_no_double_bonds = InChiParser.removeDoubleBondsSublayer(inchi2)
        return inchi1_no_double_bonds == inchi2_no_double_bonds

    def areEqualNoIsotopes(inchi1: str, inchi2: str) -> bool:
        inchi1_isotopes = InChiParser.removeIsotopicLayers(inchi1)
        inchi2_isotopes = InChiParser.removeIsotopicLayers(inchi2)
        return inchi1_isotopes == inchi2_isotopes
    
    def areEqualTautomers(inchi1: str, inchi2: str) -> bool:
        sig1 = InChiParser.getTautomerLayer(inchi1)
        sig2 = InChiParser.getTautomerLayer(inchi2)
        if sig1 is None or sig2 is None:
            return False
        return sig1 == sig2

    @staticmethod
    def get_ids(inchi1: str, inchi2: str) -> dict:
        #dict<InchiLayers, bool>: for every identity rule, returns whether true/false for each layer
            results = {}
            results[InchiLayers.COMPLETE_IDENTITY] = (
                InChi.isCompleteIdentity(inchi1, inchi2)
            )
            results[InchiLayers.INDEPENDENT_SALTS] = (
                InChi.areEqualDisolvedSalts(inchi1, inchi2) #check when removing salts: stereo
            )
            results[InchiLayers.INDEPENDENT_CHARGES] = (
                InChi.areEqualNoCharges(inchi1, inchi2)
            )
            results[InchiLayers.INDEPENDENT_DOUBLE_BONDS] = (
                InChi.areEqualNoPositionDoubleBond(inchi1, inchi2)
            )
            results[InchiLayers.TAUTOMERIC] = (
                InChi.areEqualTautomers(inchi1, inchi2)
            )
            results[InchiLayers.STEREOCHEMICAL] = (
                InChi.areEqualNoStereo(inchi1, inchi2)
            )
            results[InchiLayers.ISOTOPIC] = (
                InChi.areEqualNoIsotopes(inchi1, inchi2)
            )

            return results