import json
from pathlib import Path
from backend.inchi.config_loader import load_config
from backend.inchi.determine_levels_id import InChi

def compare_pair(inchi1, inchi2, config):
    comparison = InChi.get_ids(inchi1, inchi2, config)

    return {
        "inchi_1": inchi1,
        "inchi_2": inchi2,
        "results": {k.name: v for k, v in comparison.items()}
    }

def compare_text_files(list1, list2, config):
    results = []

    for i1 in list1:
        for i2 in list2:
            comparison = InChi.get_ids(i1, i2, config)

            results.append({
                "inchi_1": i1,
                "inchi_2": i2,
                "results": {k.name: v for k, v in comparison.items()}
            })

    return {"comparisons": results}