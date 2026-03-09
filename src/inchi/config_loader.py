from pathlib import Path
import json

DEFAULT_CONFIG_PATH = Path(__file__).resolve().parents[2] / "configs/default_config.json"

def load_config(config_path=None):

    if config_path is None:
        config_path = DEFAULT_CONFIG_PATH

    with open(config_path) as f:
        return json.load(f)