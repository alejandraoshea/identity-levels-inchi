import sys
from rdkit import Chem
from rdkit.Chem import MolToSmiles
from rdkit.Chem import inchi as rdInchi
from collections import Counter
from inchi_identity.lipid.lipid_tail_extraction import TailExtractor
from inchi_identity.lipid.lipid_pattern_generator import build_combined_patterns
from inchi_identity.inchi.inchi_parser import InChIParser
from typing import Dict, List, Optional, Tuple
from dataclasses import dataclass


@dataclass
class HeadgroupPattern:
    """
    Represents a lipid headgroup with positional specificity.
    Attributes:
        name:         Human-readable name
        smarts:       SMARTS pattern defining the headgroup
        lipid_class:  Classification (e.g. "Glycerophospholipids")
        fa_positions: sn positions where FAs attach (e.g. ["sn-1", "sn-2"])
        description:  Additional details
    """
    name: str
    smarts: str
    lipid_class: str
    fa_positions: List[str]
    description: str = ""


class LipidHeadValidator:

    HEADGROUP_PATTERNS = build_combined_patterns(
        manual_patterns={}, excel_path="Naming_Example.xlsx"
    )

    PATTERN_PRIORITY = []

    def __init__(self):
        self.compiled_patterns: Dict[str, Tuple[HeadgroupPattern, Chem.Mol]] = {}
        self.compile_patterns()

    def compile_patterns(self):
        for pattern_id, pattern in self.HEADGROUP_PATTERNS.items():
            try:
                mol_pattern = Chem.MolFromSmarts(pattern.smarts)
                if mol_pattern is not None:
                    self.compiled_patterns[pattern_id] = (pattern, mol_pattern)
                else:
                    print(f"Warning: Failed to compile pattern '{pattern_id}': {pattern.name}",
                          file=sys.stderr)
            except Exception as e:
                print(f"Error compiling pattern '{pattern_id}': {e}", file=sys.stderr)

    def matches_pattern(self, mol: Chem.Mol, pattern_id: str) -> bool:
        if pattern_id not in self.compiled_patterns:
            return False

        pattern_info, mol_pattern = self.compiled_patterns[pattern_id]

        # Step 1: SMARTS substructure match
        # Excel patterns use explicit ([H]) notation, so add Hs before matching.
        # Chirality is NOT checked here — stereo is validated separately via
        # InChI layer comparison; the SMARTS check is topology-only.
        mol_h = Chem.AddHs(mol)
        if not mol_h.HasSubstructMatch(mol_pattern, useChirality=False):
            return False

        # Step 2: FA chain guard — require at least a 2-carbon substituent so
        # pure headgroup fragments don't match, but allow short-chain test lipids.
        if pattern_info.fa_positions:
            if not LipidAnalysis.has_long_carbon_chain(mol, min_len=2):
                return False

        # Step 3: Positional constraints (pattern-specific)
        if pattern_id == "saccharolipid_o_acyl":
            ring_info = mol.GetRingInfo()
            # Build atom → rings map for fast O-ring lookup
            ring_map: Dict[int, List[set]] = {}
            for ring in ring_info.AtomRings():
                for idx in ring:
                    ring_map.setdefault(idx, []).append(set(ring))

            def in_O_ring(idx: int) -> bool:
                return any(
                    any(mol.GetAtomWithIdx(i).GetSymbol() == 'O' for i in r)
                    for r in ring_map.get(idx, [])
                )

            p_ring = Chem.MolFromSmarts("[C;R][OX2][CX3](=[OX1])[#6]~[#6]~[#6]~[#6]")
            p_exo  = Chem.MolFromSmarts("[CH2X4][OX2][CX3](=[OX1])[#6]~[#6]~[#6]~[#6]")

            o_ring_ok = False
            if p_ring:
                for match in mol.GetSubstructMatches(p_ring):
                    if in_O_ring(match[0]):
                        o_ring_ok = True
                        break
            if not o_ring_ok and p_exo:
                for match in mol.GetSubstructMatches(p_exo):
                    if any(in_O_ring(n.GetIdx())
                           for n in mol.GetAtomWithIdx(match[0]).GetNeighbors()):
                        o_ring_ok = True
                        break
            if not o_ring_ok:
                return False

            p_sph = Chem.MolFromSmarts("[NX3H2][CX4H1][CX4][OX2H1]")
            if p_sph and mol.HasSubstructMatch(p_sph):
                return False

        if pattern_id == "sterol_ester_generic":
            if mol.GetRingInfo().NumRings() < 3:
                return False
            p_phos = Chem.MolFromSmarts("[PX4]")
            if p_phos and mol.HasSubstructMatch(p_phos):
                return False

        return True

    @staticmethod
    def get_inchi(mol: Chem.Mol) -> str:
        return rdInchi.MolToInchi(mol) or ''

    @staticmethod
    def classify_with_classyfire(inchi: str):
        try:
            import requests
            inchikey = rdInchi.InchiToInchiKey(inchi)
            if not inchikey:
                return None
            resp = requests.get(
                f"https://classyfire.wishartlab.com/entities/{inchikey}.json",
                timeout=5,
                headers={"Accept": "application/json"},
            )
            if resp.status_code == 200:
                superclass = (resp.json().get("superclass") or {}).get("name", "")
                return "Lipids" in superclass
            return None
        except Exception:
            return None

    @staticmethod
    def is_lipid(inchi: str, mol) -> bool:
        cf = LipidHeadValidator.classify_with_classyfire(inchi)
        if cf is not None:
            return cf
        return LipidHeadValidator().matches_any_valid_head(mol)


    def matches_any_valid_head(
        self,
        mol: Chem.Mol,
        reference_stereo_inchi: str = None,
        reference_inchi: str = None,
    ) -> bool:
        """
        Three modes depending on which reference is supplied:

        Headgroup class only (no reference):
            Pure SMARTS topology check. Lipid class recognised; stereo and FA
            positions are ignored.

        Level B/C — chain identity (reference_inchi):
            Headgroup class + connectivity match without stereo.  Stereo layers
            (/b /t /m /s) are stripped before comparing, so only chain length,
            double bond count/position, and attachment point matter.  A molecule
            with the correct chain at the wrong sn-position fails; a molecule
            with the correct chain but inverted stereocenter passes.

        Level A — stereospecific identity (reference_stereo_inchi):
            Headgroup class + full stereo match: /b (E/Z double bonds) and
            /t /m /s (tetrahedral centres).  Requires get_stereo_layer to
            have extracted those layers from the reference molecule.
        """
        class_match = False
        for priority_group in self.PATTERN_PRIORITY:
            for pattern_id in priority_group:
                if self.matches_pattern(mol, pattern_id):
                    class_match = True
                    break
            if class_match:
                break

        if not class_match:
            manual_ids = {pid for group in self.PATTERN_PRIORITY for pid in group}
            for pattern_id in self.compiled_patterns:
                if pattern_id not in manual_ids:
                    if self.matches_pattern(mol, pattern_id):
                        class_match = True
                        break

        if not class_match:
            return False

        if reference_inchi:
            # Level B/C: compare connectivity without stereo layers so that
            # chain position / composition differences are caught but stereo
            # variants of the same chain are not rejected at this level.
            mol_inchi = LipidHeadValidator.get_inchi(mol)
            if (InChIParser.removeStereoLayers(mol_inchi)
                    != InChIParser.removeStereoLayers(reference_inchi)):
                return False

        if reference_stereo_inchi:
            # Level A: compare all stereo layers (/b /t /m /s).
            mol_stereo = InChIParser.get_stereo_layer(mol)
            if mol_stereo != reference_stereo_inchi:
                return False

        return True

    def get_matching_patterns(self, mol: Chem.Mol) -> List[HeadgroupPattern]:
        for priority_group in self.PATTERN_PRIORITY:
            for pattern_id in priority_group:
                if self.matches_pattern(mol, pattern_id):
                    pattern_info, _ = self.compiled_patterns[pattern_id]
                    return [pattern_info]
        manual_ids = {pid for group in self.PATTERN_PRIORITY for pid in group}
        for pattern_id in self.compiled_patterns:
            if pattern_id not in manual_ids:
                if self.matches_pattern(mol, pattern_id):
                    pattern_info, _ = self.compiled_patterns[pattern_id]
                    return [pattern_info]
        return []

    def identify_lipid_class(self, mol: Chem.Mol) -> List[str]:
        """Return the lipid class of the best-matching pattern, or []."""
        matches = self.get_matching_patterns(mol)
        return sorted({m.lipid_class for m in matches}) if matches else []

    def validate_structure(self, mol: Chem.Mol, verbose: bool = False) -> Dict:
        """
        Full validation of a lipid structure.
        Returns:
            {
                "is_valid":        bool,
                "lipid_class":     str or None,
                "matched_patterns":List[HeadgroupPattern],
                "fa_positions":    List[str]
            }
        """
        matches = self.get_matching_patterns(mol)
        result = {
            "is_valid":         len(matches) > 0,
            "lipid_class":      matches[0].lipid_class if matches else None,
            "matched_patterns": matches,
            "fa_positions":     matches[0].fa_positions if matches else [],
        }
        if verbose and matches:
            print(f"Valid lipid: {result['lipid_class']} / {matches[0].name}", file=sys.stderr)
        return result

    @staticmethod
    def is_valid_lipid_structure(mol: Chem.Mol) -> bool:
        return LipidHeadValidator().matches_any_valid_head(mol)

class LipidAnalysis:
    MIN_TAIL_CARBONS = 2

    @staticmethod
    def parse_smiles(smiles: str) -> Chem.Mol:
        #Parse a SMILES string into an RDKit molecule.
        if not smiles:
            return None
        if "[G]" in smiles or "[R]" in smiles:
            return None
        return Chem.MolFromSmiles(smiles.strip())

    @staticmethod
    def has_long_carbon_chain(mol, min_len=8):
        for atom in mol.GetAtoms():
            if atom.GetAtomicNum() != 6:
                continue
            visited: set = set()
            stack = [(atom.GetIdx(), 0)]
            while stack:
                idx, length = stack.pop()
                if idx in visited:
                    continue
                visited.add(idx)
                a = mol.GetAtomWithIdx(idx)
                if a.GetAtomicNum() != 6:
                    continue
                length += 1
                if length >= min_len:
                    return True
                for nbr in a.GetNeighbors():
                    stack.append((nbr.GetIdx(), length))
        return False

    @staticmethod
    def tail_sig_levelB(t):
        return (t["C"], t["DB"], t["O"], t["DB_positions"], t["O_positions"])

    @staticmethod
    def tail_sig_levelC(t):
        return (t["C"], t["DB"], t["O"])

    @staticmethod
    def atom_count(mol):
        counts: Counter = Counter()
        for atom in mol.GetAtoms():
            counts[atom.GetSymbol()] += 1
        return counts


class LipidComparator:
    @staticmethod
    def lipid_signature(mol):
        head_atoms = TailExtractor.detect_head_atoms(mol)
        if not head_atoms:
            return None
        tails = TailExtractor.extract_tails(mol)
        return tuple(sorted([
            (tail["C"], tail["DB"], tuple(sorted(tail["O_positions"])))
            for tail in tails
        ]))