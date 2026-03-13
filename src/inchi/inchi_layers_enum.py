from enum import Enum

class InchiLayers(Enum):
    COMPLETE_IDENTITY = "complete_identity"
    ISOTOPIC = "isotopic"
    INDEPENDENT_SALTS = "independent_salts"
    INDEPENDENT_CHARGES = "independent_charges"
    TAUTOMERIC = "tautomeric"
    TAUTOMERIC_1 = "tautomeric_1"
    TAUTOMERIC_2 = "tautomeric_2"
    UNKNOWN_TAUTOMERY = "unknown_tautomery"
    STEREOCHEMICAL = "stereochemical" 
    STEREOCHEMICAL_CIS_TRANS = "stereochemical_cis_trans" 
    INDEPENDENT_DOUBLE_BONDS = "independent_double_bonds"
    SUBSTITUENT_POSITION_INDEPENDENT = "substituent_position_independent"