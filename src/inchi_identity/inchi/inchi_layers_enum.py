from enum import Enum

class InchiLayers(Enum):
    COMPLETE_IDENTITY = "complete_identity"
    ISOTOPIC_INDEPENDENCE = "isotopic_independent"
    SALTS_INDEPENDENCE = "salt_independent"
    CHARGES_INDEPENDENCE = "charge_independent"
    DOUBLE_BONDS_INDEPENDENCE = "double_bond_position_independent"
    CIS_TRANS_INDEPENDENCE = "cis_trans_independent"
    SN_POSITION_INDEPENDENCE = "sn_position_independent"
    CHAIN_POSITION_INDEPENDENCE = "chain_position_independent"
    SUM_COMPOSITION_INDEPENDENCE = "sum_composition_independent"
    TAUTOMER_INDEPENDENCE = "tautomer_independent"
