class MgfParser:

    def parse_mgf(file_path):
        entries = []
        current = {}

        with open(file_path, "r") as f:
            for line in f:
                line = line.strip()

                if line == "BEGIN IONS":
                    current = {}

                elif line == "END IONS":
                    if current:
                        entries.append(current)

                elif "=" in line:
                    key, value = line.split("=", 1)
                    current[key.upper()] = value

        return entries

    def extract_inchis(entries):
        inchis = []

        for entry in entries:
            if "INCHI" in entry:
                inchis.append(entry["INCHI"])
            elif "SMILES" in entry:
                inchis.append(entry["SMILES"]) 

        return inchis
    
    def extract_structures(entries):
        result = []

        for entry in entries:
            structure = None

            if "INCHI" in entry:
                structure = entry["INCHI"]
            elif "SMILES" in entry:
                structure = entry["SMILES"]

            if structure:
                result.append({
                    "structure": structure,
                    "entry": entry
                })

        return result

"""
Simplified MGF deduplication using string-based InChI comparison.
Bypasses RDKit entirely for COMPLETE_IDENTITY level.
"""

import json
from pathlib import Path
from typing import List, Dict, Optional, Tuple
from dataclasses import dataclass, asdict
from datetime import datetime


@dataclass
class UnificationChange:
    """Tracks InChI unification changes"""
    original_inchi: str
    canonical_inchi: str
    source_file: str
    entry_title: str = ""
    timestamp: str = ""
    
    def __post_init__(self):
        if not self.timestamp:
            self.timestamp = datetime.now().isoformat()


class SimpleMgfDeduplicator:
    """
    Simple MGF deduplication using direct InChI string comparison.
    Works for COMPLETE_IDENTITY without needing RDKit.
    """
    
    def __init__(self, level: str = "COMPLETE_IDENTITY"):
        self.level = level
        self.changes_log: List[UnificationChange] = []
    
    def parse_mgf(self, file_path: str) -> List[Dict]:
        """Parse MGF file into list of entry dictionaries"""
        entries = []
        current = {}

        with open(file_path, "r") as f:
            for line in f:
                line = line.strip()

                if line == "BEGIN IONS":
                    current = {}

                elif line == "END IONS":
                    if current:
                        entries.append(current)

                elif "=" in line:
                    key, value = line.split("=", 1)
                    current[key.upper()] = value

        return entries
    
    def extract_inchi(self, entry: Dict) -> Optional[str]:
        """Extract InChI from MGF entry"""
        return entry.get("INCHI") or entry.get("SMILES")
    
    def inchis_match(self, inchi1: str, inchi2: str) -> bool:
        """
        Compare two InChIs based on the selected level.
        For COMPLETE_IDENTITY: exact string match
        """
        if not inchi1 or not inchi2:
            return False
        
        # Normalize
        inchi1 = inchi1.strip()
        inchi2 = inchi2.strip()
        
        if self.level == "COMPLETE_IDENTITY":
            return inchi1 == inchi2
        
        # Add other levels here if needed
        # For now, fallback to exact match
        return inchi1 == inchi2
    
    def deduplicate_file(
        self, 
        entries: List[Dict], 
        source_file: str
    ) -> Tuple[List[Dict], List[UnificationChange]]:
        """
        Deduplicate entries within a single file.
        
        Returns:
            (deduplicated_entries, changes_log)
        """
        groups = []
        local_changes = []
        
        for entry in entries:
            inchi = self.extract_inchi(entry)
            
            if not inchi:
                # No InChI - keep as unique
                groups.append({
                    "canonical_inchi": None,
                    "entries": [entry]
                })
                continue
            
            placed = False
            
            # Try to match with existing groups
            for group in groups:
                if group["canonical_inchi"] is None:
                    continue
                
                if self.inchis_match(inchi, group["canonical_inchi"]):
                    # Match found - add to group
                    group["entries"].append(entry)
                    
                    # Log if InChI is different (shouldn't happen for COMPLETE_IDENTITY)
                    local_changes.append(UnificationChange(
                            original_inchi=inchi,
                            canonical_inchi=group["canonical_inchi"],
                            source_file=source_file,
                            entry_title=entry.get("TITLE", entry.get("NAME", ""))
                    ))
                    
                    placed = True
                    break
            
            if not placed:
                # Create new group with this InChI as canonical
                groups.append({
                    "canonical_inchi": inchi,
                    "entries": [entry]
                })
        
        # Merge entries within each group
        deduplicated = []
        for group in groups:
            merged_entry = self.merge_entries(group["entries"])
            
            # Update to canonical InChI
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
    ) -> Tuple[List[Dict], List[UnificationChange]]:
        """
        M×N comparison between two files.
        Use entries from A as canonical.
        """
        unified = list(entries_a)  # Start with all from A
        cross_changes = []
        
        for entry_b in entries_b:
            inchi_b = self.extract_inchi(entry_b)
            
            if not inchi_b:
                # No InChI - add as unique
                unified.append(entry_b)
                continue
            
            matched = False
            
            # Compare against all entries in A
            for i, entry_a in enumerate(unified):
                inchi_a = self.extract_inchi(entry_a)
                
                if not inchi_a:
                    continue
                
                if self.inchis_match(inchi_b, inchi_a):
                    merged = self.merge_entries([entry_a, entry_b])
                    
                    merged["INCHI"] = inchi_a
                    
                    cross_changes.append(UnificationChange(
                            original_inchi=inchi_b,
                            canonical_inchi=inchi_a,
                            source_file=f"{source_b} → {source_a}",
                            entry_title=entry_b.get("TITLE", entry_b.get("NAME", ""))
                    ))
                    
                    unified[i] = merged
                    matched = True
                    break
            
            if not matched:
                unified.append(entry_b)
        
        return unified, cross_changes
    
    def merge_entries(self, entries: List[Dict]) -> Dict:
        """
        Merge multiple MGF entries.
        Combines metadata; first entry's values take precedence.
        """
        if len(entries) == 1:
            return entries[0].copy()
        
        merged = {}
        
        all_keys = set()
        for entry in entries:
            all_keys.update(entry.keys())
        
        for key in all_keys:
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
    
    def write_mgf(self, entries: List[Dict], output_path: str):
        """Write entries to MGF file"""
        with open(output_path, "w") as f:
            for entry in entries:
                f.write("BEGIN IONS\n")
                for key, value in entry.items():
                    f.write(f"{key}={value}\n")
                f.write("END IONS\n\n")
    
    def write_log(self, output_path: str):
        """Write unification log to file"""
        log_data = {
            "level": self.level,
            "timestamp": datetime.now().isoformat(),
            "total_changes": len(self.changes_log),
            "changes": [asdict(change) for change in self.changes_log]
        }
        
        with open(output_path, "w") as f:
            json.dump(log_data, f, indent=2)
    
    def process_files(
        self,
        file_path_a: str,
        file_path_b: str,
        output_mgf: Optional[str] = None,
        output_log: Optional[str] = None
    ) -> Dict:
        """
        Complete workflow: deduplicate → cross-compare → merge → output.
        """
        
        entries_a = self.parse_mgf(file_path_a)
        entries_b = self.parse_mgf(file_path_b)
        
        source_a = Path(file_path_a).name
        source_b = Path(file_path_b).name
        
        # Step 1: Deduplicate within File A
        dedup_a, changes_a = self.deduplicate_file(entries_a, source_a)
        self.changes_log.extend(changes_a)
        
        # Step 2: Deduplicate within File B
        dedup_b, changes_b = self.deduplicate_file(entries_b, source_b)
        self.changes_log.extend(changes_b)
        
        # Step 3: Cross-deduplicate (M×N)
        unified, cross_changes = self.cross_deduplicate(
            dedup_a, dedup_b, source_a, source_b
        )
        self.changes_log.extend(cross_changes)
        
        # Step 4: Write outputs
        if output_mgf:
            self.write_mgf(unified, output_mgf)
        
        if output_log:
            self.write_log(output_log)
        
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
            "level": self.level
        }