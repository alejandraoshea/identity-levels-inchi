from enum import Enum

class InchiLayers(Enum):
    COMPLETE_IDENTITY = "complete_identity"
    INDEPENDENT_SALTS = "independent_salts"
    INDEPENDENT_CHARGES = "independent_charges"
    INDEPENDENT_DOUBLE_BONDS = "independent_double_bonds"
    TAUTOMERIC = "tautomeric"
    STEREOCHEMICAL = "stereochemical" #extended SMILES
    ISOTOPIC = "isotopic" #presencia heavy metals
    
