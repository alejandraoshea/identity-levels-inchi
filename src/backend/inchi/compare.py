import json, tempfile, os
from pathlib import Path
from backend.inchi.determine_levels_id import InChI
from backend.parsers.mgf_parser import MgfParser, SimpleMgfDeduplicator
from typing import List, Dict, Optional, Tuple
from dataclasses import dataclass, asdict
from datetime import datetime


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
            if res is not None:
                results.append(res)

    elif mode == "cross":
        for i1 in list1:
            for i2 in list2:
                res = process(i1, i2)
                if res is not None:
                    results.append(res)

    return {"comparisons": results}

@dataclass
class IdentityChange:
    """Tracks InChI identity transformations"""
    original_inchi: str
    canonical_inchi: str
    level: str
    source_file: str
    entry_title: str = ""
    timestamp: str = ""
    
    def __post_init__(self):
        if not self.timestamp:
            self.timestamp = datetime.now().isoformat()


class MgfDeduplicator:
    """
    Handles MGF file deduplication and merging based on identity levels.
    """
    
    def __init__(self, config: dict, level: str = "COMPLETE_IDENTITY"):
        """
        Args:
            config: Identity configuration
            level: Identity level to use for comparison
        """
        self.config = config
        self.level = level
        self.changes_log: List[IdentityChange] = []
    
    def deduplicate_file(
        self, 
        entries: List[Dict], 
        source_file: str
    ) -> Tuple[List[Dict], List[IdentityChange]]:
        """
        Deduplicate entries within a single MGF file.
        
        Args:
            entries: List of MGF entries (dicts)
            source_file: Name of source file (for logging)
        
        Returns:
            (deduplicated_entries, changes_log)
        """
        groups = []
        local_changes = []
        
        for entry in entries:
            inchi = self._extract_inchi(entry)
            if not inchi:
                groups.append({
                    "canonical_inchi": None,
                    "entries": [entry]
                })
                continue
            
            placed = False
            
            for group in groups:
                if group["canonical_inchi"] is None:
                    continue
                
                comparison = InChI.get_ids(
                    inchi,
                    group["canonical_inchi"],
                    self.config
                )
                
                if comparison.get(self.level):
                    group["entries"].append(entry)
                    
                    if inchi != group["canonical_inchi"]:
                        local_changes.append(IdentityChange(
                            original_inchi=inchi,
                            canonical_inchi=group["canonical_inchi"],
                            level=self.level,
                            source_file=source_file,
                            entry_title=entry.get("TITLE", entry.get("NAME", ""))
                        ))
                    
                    placed = True
                    break
            
            if not placed:
                groups.append({
                    "canonical_inchi": inchi,
                    "entries": [entry]
                })
        
        deduplicated = []
        for group in groups:
            merged_entry = self._merge_entries(group["entries"])
            
            if group["canonical_inchi"]:
                merged_entry["INCHI"] = group["canonical_inchi"]
            
            deduplicated.append(merged_entry)
        
        return deduplicated, local_changes
    
    def cross_deduplicate(
        self,
        entries_a: List[Dict],
        entries_b: List[Dict],
        source_a: str = "File A",
        source_b: str = "File B"
    ) -> Tuple[List[Dict], List[IdentityChange]]:
        """
        Cross-deduplicate between two lists of entries (already deduplicated internally).
        
        Args:
            entries_a: Deduplicated entries from file A
            entries_b: Deduplicated entries from file B
            source_a: Name of file A
            source_b: Name of file B
        
        Returns:
            (unified_entries, cross_changes_log)
        """
        unified = list(entries_a) 
        cross_changes = []
        
        for entry_b in entries_b:
            inchi_b = self._extract_inchi(entry_b)
            if not inchi_b:
                unified.append(entry_b)
                continue
            
            matched = False
            
            for i, entry_a in enumerate(unified):
                inchi_a = self._extract_inchi(entry_a)
                if not inchi_a:
                    continue
                
                comparison = InChI.get_ids(
                    inchi_b,
                    inchi_a,
                    self.config
                )
                
                if comparison.get(self.level):
                    merged = self._merge_entries([entry_a, entry_b])
                    
                    merged["INCHI"] = inchi_a
                    
                    if inchi_b != inchi_a:
                        cross_changes.append(IdentityChange(
                            original_inchi=inchi_b,
                            canonical_inchi=inchi_a,
                            level=self.level,
                            source_file=f"{source_b} → {source_a}",
                            entry_title=entry_b.get("TITLE", entry_b.get("NAME", ""))
                        ))
                    
                    unified[i] = merged
                    matched = True
                    break
            
            if not matched:
                unified.append(entry_b)
        
        return unified, cross_changes
    
    def process_files(
        self,
        file_path_a: str,
        file_path_b: str,
        output_mgf: Optional[str] = None,
        output_log: Optional[str] = None
    ) -> Dict:
        """
        Complete workflow: deduplicate → cross-compare → merge → output.
        
        Args:
            file_path_a: Path to first MGF file
            file_path_b: Path to second MGF file
            output_mgf: Path for output MGF file (optional)
            output_log: Path for changes log JSON (optional)
        
        Returns:
            Dictionary with results and statistics
        """
        entries_a = MgfParser.parse_mgf(file_path_a)
        entries_b = MgfParser.parse_mgf(file_path_b)
        
        source_a = Path(file_path_a).name
        source_b = Path(file_path_b).name
    
        # Step 1: Deduplicate within each file
        dedup_a, changes_a = self.deduplicate_file(entries_a, source_a)
        dedup_b, changes_b = self.deduplicate_file(entries_b, source_b)

        self.changes_log.extend(changes_a)
        self.changes_log.extend(changes_b)
        
        # Step 2: Cross-deduplicate between files
        unified, cross_changes = self.cross_deduplicate(
            dedup_a, dedup_b, source_a, source_b
        )

        self.changes_log.extend(cross_changes)
        
        # Step 3: Write output files
        if output_mgf:
            self._write_mgf(unified, output_mgf)
        
        if output_log:
            self._write_log(output_log)
        
        return {
            "input_counts": {
                source_a: len(entries_a),
                source_b: len(entries_b)
            },
            "deduplicated_counts": {
                source_a: len(dedup_a),
                source_b: len(dedup_b)
            },
            "unified_count": len(unified),
            "changes_count": len(self.changes_log),
            "changes_log": [asdict(change) for change in self.changes_log],
            "level": self.level,
            "output_mgf": output_mgf,
            "output_log": output_log
        }
    
    def _extract_inchi(self, entry: Dict) -> Optional[str]:
        """Extract InChI from MGF entry"""
        return entry.get("INCHI") or entry.get("SMILES")
    
    def _merge_entries(self, entries: List[Dict]) -> Dict:
        """
        Merge multiple MGF entries into one.
        Combines metadata and MS/MS spectra.
        """
        if len(entries) == 1:
            return entries[0].copy()
        
        merged = {}
        
        for key in entries[0].keys():
            values = []
            for entry in entries:
                if key in entry:
                    val = entry[key]
                    if val not in values:
                        values.append(val)
            
            if len(values) == 1:
                merged[key] = values[0]
            elif len(values) > 1:
                if key in ["TITLE", "NAME"]:
                    merged[key] = " | ".join(str(v) for v in values)
                else:
                    merged[key] = values[0]  
        
        
        return merged
    
    def _write_mgf(self, entries: List[Dict], output_path: str):
        """Write entries to MGF file"""
        with open(output_path, "w") as f:
            for entry in entries:
                f.write("BEGIN IONS\n")
                for key, value in entry.items():
                    f.write(f"{key}={value}\n")
                f.write("END IONS\n\n")
    
    def _write_log(self, output_path: str):
        """Write changes log to JSON"""
        log_data = {
            "level": self.level,
            "timestamp": datetime.now().isoformat(),
            "total_changes": len(self.changes_log),
            "changes": [asdict(change) for change in self.changes_log]
        }
        
        with open(output_path, "w") as f:
            json.dump(log_data, f, indent=2)


def compare_mgf_files(
    file1,
    file2,
    config,
    level="COMPLETE_IDENTITY",
    merge_msms=False,
    output_mgf=None,
    output_log=None
):
    """
    Compare and deduplicate MGF files with unified output.
    
    Args:
        file1: Path to first MGF file or list of entries
        file2: Path to second MGF file or list of entries
        config: Identity configuration
        level: Identity level for comparison (default: COMPLETE_IDENTITY)
        merge_msms: If True, use unified deduplication workflow
        output_mgf: Optional output path for unified MGF
        output_log: Optional output path for changes log
    
    Returns:
        Dictionary with results
    """
    if merge_msms:
        deduplicator = MgfDeduplicator(config, level)
        
        if isinstance(file1, str):
            return deduplicator.process_files(
                file1, file2,
                output_mgf=output_mgf,
                output_log=output_log
            )
        else:
            source_a = "File A"
            source_b = "File B"
            
            dedup_a, changes_a = deduplicator.deduplicate_file(file1, source_a)
            dedup_b, changes_b = deduplicator.deduplicate_file(file2, source_b)
            
            deduplicator.changes_log.extend(changes_a)
            deduplicator.changes_log.extend(changes_b)
            
            unified, cross_changes = deduplicator.cross_deduplicate(
                dedup_a, dedup_b, source_a, source_b
            )
            
            deduplicator.changes_log.extend(cross_changes)
            
            return {
                "input_counts": {
                    source_a: len(file1),
                    source_b: len(file2)
                },
                "deduplicated_counts": {
                    source_a: len(dedup_a),
                    source_b: len(dedup_b)
                },
                "unified_count": len(unified),
                "changes_count": len(deduplicator.changes_log),
                "changes_log": [asdict(change) for change in deduplicator.changes_log],
                "unified_entries": unified,
                "level": level
            }
    
    else:
        if isinstance(file1, str):
            entries1 = MgfParser.parse_mgf(file1)
            entries2 = MgfParser.parse_mgf(file2)
        else:
            entries1 = file1
            entries2 = file2

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

    return merged

def compare_mgf_files_simple(
    file1: str,
    file2: str,
    level: str = "COMPLETE_IDENTITY",
    output_mgf: Optional[str] = None,
    output_log: Optional[str] = None
) -> Dict:
    """
    Simple MGF comparison using string-based InChI matching.
    Bypasses RDKit completely.
    
    Args:
        file1: Path to first MGF file or list of entries
        file2: Path to second MGF file or list of entries
        level: Identity level (default: COMPLETE_IDENTITY)
        output_mgf: Output path for unified MGF
        output_log: Output path for unification log
    
    Returns:
        Dictionary with results
    """
    deduplicator = SimpleMgfDeduplicator(level=level)
    
    if isinstance(file1, str):
        return deduplicator.process_files(
            file1, file2,
            output_mgf=output_mgf,
            output_log=output_log
        )
    else:
        with tempfile.NamedTemporaryFile(mode='w', delete=False, suffix='.mgf') as tmp1:
            deduplicator.write_mgf(file1, tmp1.name)
            tmp1_path = tmp1.name
        
        with tempfile.NamedTemporaryFile(mode='w', delete=False, suffix='.mgf') as tmp2:
            deduplicator.write_mgf(file2, tmp2.name)
            tmp2_path = tmp2.name
        
        result = deduplicator.process_files(
            tmp1_path, tmp2_path,
            output_mgf=output_mgf,
            output_log=output_log
        )
        
        os.unlink(tmp1_path)
        os.unlink(tmp2_path)
        
        return result