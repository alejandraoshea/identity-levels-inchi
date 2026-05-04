from backend.input.input_inchi import InputInChI
from backend.inchi.determine_levels_id import InChI
from backend.inchi.config_loader import load_config
from backend.inchi import InChI
from backend.input.input_inchi import InputInChI

config = load_config()

InputInChI.input_inchi(InChI.get_ids, config)
InputInChI.input_inchi(InChI.get_ids)
