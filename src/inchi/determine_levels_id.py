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
        return(InChiParser.getMainLayer(inchi1) == InChiParser.getMainLayer(inchi2) and
            InChiParser.getAtomConnectionsSublayer(inchi1) == InChiParser.getAtomConnectionsSublayer(inchi2) and
            InChiParser.getHydrogenAtomsSublayer(inchi1) == InChiParser.getHydrogenAtomsSublayer(inchi2)
        )
    
    #charge layer (sublayers: charge and proton)
    def charge_layer(inchi1: str, inchi2: str) -> bool:
        return (
            InChiParser.getChargeSublayer(inchi1) == InChiParser.getChargeSublayer(inchi2) and
            InChiParser.getProtonSublayer(inchi1) == InChiParser.getProtonSublayer(inchi2)
        )
    
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
        except Exception:
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

    def areEqualNoCharges(inchi1: str, inchi2:str) -> bool:
        inchi1_no_charge = InChiParser.removeChargeLayersUsingParser(inchi1)
        inchi2_no_charge = InChiParser.removeChargeLayersUsingParser(inchi2)
        return inchi1_no_charge == inchi2_no_charge

    def areEqualNoPositionDoubleBond(inchi1: str, inchi2:str) -> bool:
        mol1 = InChi.mol_from_inchi(inchi1)
        mol2 = InChi.mol_from_inchi(inchi2)
        if not mol1 or not mol2:
            return False

        mcs = rdFMCS.FindMCS(
            [mol1, mol2],
            atomCompare=rdFMCS.AtomCompare.CompareElements,
            bondCompare=rdFMCS.BondCompare.CompareOrder
        )
        return mcs.numAtoms == mol1.GetNumAtoms() == mol2.GetNumAtoms()


    def areEqualNoIsotopes(inchi1: str, inchi2: str) -> bool:
        mol1 = InChi.mol_from_inchi(inchi1)
        mol2 = InChi.mol_from_inchi(inchi2)
        
        if not mol1 or not mol2:
            return False
        
        mcs = rdFMCS.FindMCS(
            [mol1, mol2],
            atomCompare=rdFMCS.AtomCompare.CompareElements,  
            bondCompare=rdFMCS.BondCompare.CompareOrderExact
        )
        
        return mcs.numAtoms == mol1.GetNumAtoms() == mol2.GetNumAtoms()

    def get_ids(inchi1: str, inchi2: str) -> dict:
            """
            Returns: dict<InchiLayers, bool>
            For every identity rule, returns whether it is satisfied by the pair.
            """

            results = {}

            results[InchiLayers.COMPLETE_IDENTITY] = (
                InChi.isCompleteIdentity(inchi1, inchi2)
            )

            results[InchiLayers.INDEPENDENT_SALTS] = (
                InChi.equal_ignoring_salts(inchi1, inchi2)
            )

            results[InchiLayers.INDEPENDENT_CHARGES] = (
                InChi.main_layer(inchi1, inchi2) and
                not InChi.charge_layer(inchi1, inchi2)
            )

            results[InchiLayers.INDEPENDENT_DOUBLE_BONDS] = (
                InChi.equal_ignoring_double_bond_position(inchi1, inchi2)
            )

            results[InchiLayers.TAUTOMERIC] = (
                InChi.main_layer(inchi1, inchi2) and
                InChiParser.getHydrogenAtomsSublayer(inchi1)
                != InChiParser.getHydrogenAtomsSublayer(inchi2)
            )

            results[InchiLayers.STEREOCHEMICAL] = (
                InChi.main_layer(inchi1, inchi2) and
                not InChi.stereochemical_layer(inchi1, inchi2)
            )

            results[InchiLayers.ISOTOPIC] = (
                InChi.equal_ignoring_isotopes(inchi1, inchi2)
            )

            return results


#TODO:
#MAP<IDS,BOOLEAN>GETIDS(INCHI1,INCHI2)
