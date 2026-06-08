from inchi_identity.input.input_inchi import InputInChI
from inchi_identity.inchi.determine_levels_id import InChI
from inchi_identity.inchi.config_loader import load_config
from inchi_identity.inchi import InChI
from inchi_identity.input.input_inchi import InputInChI

config = load_config()

InputInChI.input_inchi(InChI.get_ids, config)
InputInChI.input_inchi(InChI.get_ids)
