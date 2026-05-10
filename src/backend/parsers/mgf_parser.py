import json
from pathlib import Path
from typing import List, Dict, Optional
from dataclasses import dataclass, asdict

from backend.inchi.determine_levels_id import InChI
from backend.inchi.inchi_layers_enum import InchiLayers

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
    

@dataclass
class UnificationChange:
    """Tracks InChI unification changes"""
    original_inchi: str
    canonical_inchi: str

class SimpleMgfDeduplicator:
    def __init__(self, level: str = "COMPLETE_IDENTITY", config: dict = None):
        self.level = level
        self.config = config
        self.changes_log: List[UnificationChange] = []
        
        self.level_map = {
            "COMPLETE_IDENTITY": InchiLayers.COMPLETE_IDENTITY,
            "ISOTOPIC_INDEPENDENCE": InchiLayers.ISOTOPIC_INDEPENDENCE,
            "SALTS_INDEPENDENCE": InchiLayers.SALTS_INDEPENDENCE,
            "CHARGES_INDEPENDENCE": InchiLayers.CHARGES_INDEPENDENCE,
            "DOUBLE_BONDS_INDEPENDENCE": InchiLayers.DOUBLE_BONDS_INDEPENDENCE,
            "STEREOCHEMICAL_CIS_TRANS_INDEPENDENCE": InchiLayers.STEREOCHEMICAL_CIS_TRANS_INDEPENDENCE,
            "TAUTOMER_INDEPENDENCE": InchiLayers.TAUTOMER_INDEPENDENCE,
            "SUBSTITUENT_POSITION_INDEPENDENCE": InchiLayers.SUBSTITUENT_POSITION_INDEPENDENCE
        }
    
    def parse_mgf(self, file_path: str) -> List[Dict]:
        entries = []
        current = {}
        current_peaks = []

        with open(file_path, "r") as f:
            for line in f:
                line_stripped = line.strip()

                if line_stripped == "BEGIN IONS":
                    current = {}
                    current_peaks = []

                elif line_stripped == "END IONS":
                    if current:
                        if current_peaks:
                            current["_PEAKS"] = current_peaks
                        entries.append(current)

                elif "=" in line_stripped:
                    key, value = line_stripped.split("=", 1)
                    current[key.upper()] = value
                
                else:
                    if line_stripped and current:
                        current_peaks.append(line_stripped)

        return entries
    
    def extract_inchi(self, entry: Dict) -> Optional[str]:
        return entry.get("INCHI") or entry.get("SMILES")
    
    def inchis_match(self, inchi1: str, inchi2: str) -> bool:
        if not inchi1 or not inchi2:
            return False
        
        inchi1 = inchi1.strip()
        inchi2 = inchi2.strip()
        
        if self.level == "COMPLETE_IDENTITY":
            return inchi1 == inchi2
        
        if not self.config:
            print(f"Warning: Config required for level {self.level}, falling back to COMPLETE_IDENTITY")
            return inchi1 == inchi2
        
        try:
            comparison = InChI.get_ids(inchi1, inchi2, self.config)
            level_enum = self.level_map.get(self.level, InchiLayers.COMPLETE_IDENTITY)
            
            return comparison.get(level_enum, False)
        
        except Exception as e:
            print(f"Error comparing InChIs at level {self.level}: {e}")
            return inchi1 == inchi2
    
    def unify_inchis_in_file(self, entries: List[Dict], source_file: str) -> List[Dict]:
        canonical_map = {} 
        
        for entry in entries:
            inchi = self.extract_inchi(entry)
            if not inchi:
                continue
            
            found_canonical = None
            for canonical_inchi in canonical_map.values():
                if self.inchis_match(inchi, canonical_inchi):
                    found_canonical = canonical_inchi
                    break
            
            if found_canonical:
                if inchi != found_canonical:
                    self.changes_log.append(UnificationChange(
                        original_inchi=inchi,
                        canonical_inchi=found_canonical
                    ))
                
                canonical_map[inchi] = found_canonical
            else:
                canonical_map[inchi] = inchi
        
        modified_entries = []
        for entry in entries:
            entry_copy = entry.copy()
            inchi = self.extract_inchi(entry_copy)
            
            if inchi and inchi in canonical_map:
                if "INCHI" in entry_copy:
                    entry_copy["INCHI"] = canonical_map[inchi]
                elif "SMILES" in entry_copy:
                    entry_copy["SMILES"] = canonical_map[inchi]
            
            modified_entries.append(entry_copy)
        
        return modified_entries
    
    def cross_unify(
        self,
        entries_a: List[Dict],
        entries_b: List[Dict],
        source_a: str = "File A",
        source_b: str = "File B"
    ) -> List[Dict]:
        canonical_from_a = {}
        for entry in entries_a:
            inchi = self.extract_inchi(entry)
            if inchi:
                canonical_from_a[inchi] = inchi
        
        modified_b = []
        for entry in entries_b:
            entry_copy = entry.copy()
            inchi_b = self.extract_inchi(entry_copy)
            
            if inchi_b:
                matched_canonical = None
                for inchi_a in canonical_from_a.keys():
                    if self.inchis_match(inchi_b, inchi_a):
                        matched_canonical = inchi_a
                        break
                
                if matched_canonical:
                    if inchi_b != matched_canonical:
                        self.changes_log.append(UnificationChange(
                            original_inchi=inchi_b,
                            canonical_inchi=matched_canonical
                        ))
                    
                    if "INCHI" in entry_copy:
                        entry_copy["INCHI"] = matched_canonical
                    elif "SMILES" in entry_copy:
                        entry_copy["SMILES"] = matched_canonical
            
            modified_b.append(entry_copy)
        
        return entries_a + modified_b
    
    
    def write_mgf(self, entries: List[Dict], output_path: str):
        with open(output_path, "w") as f:
            for entry in entries:
                f.write("BEGIN IONS\n")
                
                for key, value in entry.items():
                    if key != "_PEAKS":
                        f.write(f"{key}={value}\n")
                
                if "_PEAKS" in entry:
                    for peak_line in entry["_PEAKS"]:
                        f.write(f"{peak_line}\n")
                
                f.write("END IONS\n\n")
    
    def write_log(self, output_path: str):
        log_data = {
            "level": self.level,
            "total_changes": len(self.changes_log),
            "changes": [asdict(change) for change in self.changes_log]
        }
        
        with open(output_path, "w") as f:
            json.dump(log_data, f, indent=2)
    
    def process_files(
        self,file_path_a: str,file_path_b: str,
        output_mgf: Optional[str] = None,output_log: Optional[str] = None) -> Dict:
        
        entries_a = self.parse_mgf(file_path_a)
        entries_b = self.parse_mgf(file_path_b)
        
        source_a = Path(file_path_a).name
        source_b = Path(file_path_b).name
        
        changes_before = len(self.changes_log)
        entries_a_unified = self.unify_inchis_in_file(entries_a, source_a)
        changes_a = len(self.changes_log) - changes_before
        
        changes_before = len(self.changes_log)
        entries_b_unified = self.unify_inchis_in_file(entries_b, source_b)
        changes_b = len(self.changes_log) - changes_before
        
        changes_before = len(self.changes_log)
        all_entries = self.cross_unify(
            entries_a_unified, entries_b_unified, source_a, source_b
        )
        changes_cross = len(self.changes_log) - changes_before
        
        if output_mgf:
            self.write_mgf(all_entries, output_mgf)

        if output_log:
            self.write_log(output_log)
            
        return {
            "input_counts": {
                source_a: len(entries_a),
                source_b: len(entries_b)
            },
            "output_count": len(all_entries),
            "changes_count": len(self.changes_log),
            "changes_breakdown": {
                f"{source_a}_internal": changes_a,
                f"{source_b}_internal": changes_b,
                f"{source_b}_to_{source_a}": changes_cross
            },
            "changes_log": [asdict(change) for change in self.changes_log],
            "level": self.level
        }
