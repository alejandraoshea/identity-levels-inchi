import json
from pathlib import Path
from backend.inchi.determine_levels_id import InChI
from backend.parsers.mgf_parser import MgfParser


def compare_pair(inchi1, inchi2, config):
    comparison = InChI.get_ids(inchi1, inchi2, config)

    return {
        "inchi_1": inchi1,
        "inchi_2": inchi2,
        "results": {k.name: bool(v) for k, v in comparison.items()}
    }


def read_file_lines(file_path):
    return [
        line.strip()
        for line in Path(file_path).read_text().splitlines()
        if line.strip()
    ]


def extract_matches(comparison):
    return {
        k.name: True
        for k, v in comparison.items()
        if bool(v)
    }


def compare_text_files(list1, list2, config, mode="pairwise", only_equal=False):
    results = []

    def process(i1, i2):
        comparison = InChI.get_ids(i1, i2, config)

        if only_equal:
            matches = extract_matches(comparison)
            if not matches:
                return None
            return {
                "inchi_1": i1,
                "inchi_2": i2,
                "matches": matches
            }
        else:
            return {
                "inchi_1": i1,
                "inchi_2": i2,
                "results": {k.name: bool(v) for k, v in comparison.items()}
            }

    if mode == "pairwise":
        for i1, i2 in zip(list1, list2):
            res = process(i1, i2)
            if res is not None:        # ← use "is not None", not truthiness
                results.append(res)

    elif mode == "cross":
        for i1 in list1:
            for i2 in list2:
                res = process(i1, i2)
                if res is not None:
                    results.append(res)

    return {"comparisons": results}
    
def compare_mgf_files(
    file1,
    file2,
    config,
    level="COMPLETE_IDENTITY",
    merge_msms=False
):
    entries1 = MgfParser.parse_mgf(file1)
    entries2 = MgfParser.parse_mgf(file2)

    structs1 = MgfParser.extract_structures(entries1)
    structs2 = MgfParser.extract_structures(entries2)

    all_structs = structs1 + structs2

    groups = []

    for item in all_structs:
        placed = False

        for group in groups:
            rep = group["representative"]

            comparison = InChI.get_ids(
                item["structure"],
                rep,
                config
            )

            if comparison.get(level):
                group["entries"].append(item["entry"])
                placed = True
                break

        if not placed:
            groups.append({
                "representative": item["structure"],
                "entries": [item["entry"]]
            })

    if merge_msms:
        unified = []

        for group in groups:
            merged_entry = merge_group_entries(group["entries"])
            unified.append({
                "representative": group["representative"],
                "merged_entry": merged_entry,
                "num_spectra": len(group["entries"])
            })

        return {"groups": unified}

    return {"groups": groups}

def merge_group_entries(entries):
    merged = {}
    base = entries[0].copy()

    for key in base:
        values = set()

        for e in entries:
            if key in e:
                values.add(e[key])

        merged[key] = list(values) if len(values) > 1 else list(values)[0]

    # DEFINIR para añadir (eg:
    # - merge de picos MS/MS
    # - normalización
    # - clustering)

    return merged