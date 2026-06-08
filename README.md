# inchi-identity

A hierarchical molecular identity comparison framework for untargeted metabolomics,
implemented as a Python CLI and published on PyPI.

Instead of binary exact matching, this tool evaluates molecular equivalence across six
progressive normalization layers, returning a structured equivalence profile that reflects
the structural resolution actually supported by the experimental evidence. For lipids, a
dedicated four-level hierarchy (Levels A–D) addresses cis/trans geometry, sn-position,
intra-chain double bond position, and global sum composition.

Both InChI strings and SMILES are accepted as input — SMILES are automatically converted
to canonical InChI before comparison.

---

## Comparison Layers

| Layer | Name | Description |
|-------|------|-------------|
| 1 | Complete Identity | Exact InChI string equality |
| 2 | Isotopic Independence | Equality after `/i` layer removal |
| 3 | Salt Independence | Equality after counterion removal (RDKit SaltRemover) |
| 4 | Charge Independence | Equality after charge/protonation normalization |
| 5 | Stereochemical/Isomeric Independence | Levels A–D (see below) |
| 6 | Tautomeric Independence | Equality after canonical tautomer generation |

Preprocessing is **cumulative**: each layer applies all normalizations from preceding
layers before performing its own comparison.

### Lipid hierarchy (Layer 5, Levels A–D)

For molecules classified as lipids, Layer 5 applies a four-level structural abstraction
cascade. Level A applies to all molecule types; Levels B–D apply to lipids only.

| Level | Name | Description |
|-------|------|-------------|
| A | Cis/Trans Independence | Removes `/b` double-bond geometry layer |
| B | sn-Position Independence | Compares acyl chains as an unordered set |
| C | Intra-Chain Position Independence | Discards double bond and substituent positions within each chain |
| D | Global Composition | Retains only total carbon count and double bond count across the whole molecule |

Lipid classification uses the ClassyFire API, with an RDKit-based heuristic fallback.
Headgroup validation uses a library of 295 SMARTS patterns (40 manually curated + 255
generated at runtime by combining 15 glycolipid templates with 17 monosaccharide patterns).

---

## Installation

RDKit must be installed via conda before installing this package:

```bash
git clone https://github.com/alejandraoshea/identity-levels-inchi.git
cd identity-levels-inchi
conda env create -f conda_env.yml
conda activate inchi-identity
pip install -e .
```

Or install directly from PyPI:

```bash
pip install inchi-identity
```

Verify the installation:

```bash
inchi --help
```

### InChI Trust executable (optional, required for Layer 6)

Layer 6 uses the InChI Trust executable for canonical tautomer generation. Download it
from https://www.inchi-trust.org/downloads/ and set the path:

```bash
export INCHITRUST_PATH=/path/to/inchi-1
```

If the executable is not available, Layer 6 automatically falls back to RDKit's
`TautomerEnumerator`.

---

## CLI Usage

### Mode 1 — Compare two structures across all layers

```bash
inchi compare-pair "<inchi_1>" "<inchi_2>"
inchi compare-pair "<inchi_1>" "<inchi_2>" > result.json
```

SMILES are also accepted:

```bash
inchi compare-pair "CC(=O)O" "CC(=O)[O-].[Na+]"
```

### Mode 2 — Compare with selected layers only

```bash
inchi compare-pair-layers "<inchi_1>" "<inchi_2>" \
    --layers isotope charge double_bond tautomer
```

Available layer names: `complete_identity`, `isotope`, `salt`, `charge`, `double_bond`, `tautomer`

### Mode 3 — File-based comparison

Input files contain one InChI or SMILES per line. SMILES and InChI can be mixed.

**Pairwise** (entry *i* from file 1 vs entry *i* from file 2):

```bash
inchi compare file1.txt file2.txt > result_pairwise.json
```

**Cross-comparison** (all vs all, n×m pairs):

```bash
inchi compare file1.txt file2.txt --mode cross > result_cross.json
```

**Filter to equivalent pairs only** (works with both modes):

```bash
inchi compare file1.txt file2.txt --only-equal
inchi compare file1.txt file2.txt --mode cross --only-equal
```

When `--only-equal` is active, the output key changes from `results` to `matches` and
only layers that evaluated to `true` are included.

### Mode 4 — MGF spectral library unification

Normalizes molecular identifiers embedded in two MGF files at a chosen equivalence layer,
producing a unified MGF file and a JSON change log. InChI and SMILES identifiers are
detected automatically per entry.

```bash
inchi compare-mgf file1.mgf file2.mgf \
    --layer CHARGES_INDEPENDENCE \
    --output-mgf unified.mgf \
    --output-log unified_log.json
```

Available `--layer` values:

```
COMPLETE_IDENTITY
ISOTOPIC_INDEPENDENCE
SALTS_INDEPENDENCE
CHARGES_INDEPENDENCE
DOUBLE_BONDS_INDEPENDENCE
STEREOCHEMICAL_CIS_TRANS_INDEPENDENCE
TAUTOMER_INDEPENDENCE
```

The pipeline: (1) parses both files and extracts identifiers, (2) internally normalizes
each file so equivalent entries share a canonical InChI, (3) performs a cross-file
comparison and replaces identifiers in file 2 with the canonical form from file 1,
(4) writes the unified MGF and a JSON log recording every transformation. Entries without
identifiers and all spectral peak data are preserved unchanged.

---

## Output format

### Pairwise / file comparison

```json
{
    "inchi_1": "InChI=1S/C5H11NO2/c1-6(2,3)4-5(7)8/h4H2,1-3H3/p+1",
    "inchi_2": "InChI=1S/C5H11NO2/c1-6(2,3)4-5(7)8/h4H2,1-3H3",
    "results": {
        "COMPLETE_IDENTITY": false,
        "ISOTOPIC_INDEPENDENCE": false,
        "SALTS_INDEPENDENCE": false,
        "CHARGES_INDEPENDENCE": true,
        "DOUBLE_BONDS_INDEPENDENCE": false,
        "STEREOCHEMICAL_CIS_TRANS_INDEPENDENCE": false,
        "TAUTOMER_INDEPENDENCE": false
    }
}
```

### With `--only-equal`

Only layers that evaluated to `true` are shown, under the `matches` key:

```json
{
    "inchi_1": "InChI=1S/C5H11NO2/c1-6(2,3)4-5(7)8/h4H2,1-3H3/p+1",
    "inchi_2": "InChI=1S/C5H11NO2/c1-6(2,3)4-5(7)8/h4H2,1-3H3",
    "matches": {
        "CHARGES_INDEPENDENCE": true
    }
}
```

### MGF change log

```json
{
    "layer": "CHARGES_INDEPENDENCE",
    "total_changes": 3,
    "changes": [
        {
            "original_structure": "original_structure": "InChI=1S/C10H16N5O13P3/c11-8-5-9(13-2-12-8)15(3-14-5)10-7(17)6(16)4(26-10)1-25-30(21,22)28-31(23,24)27-29(18,19)20/h2-4,6-7,10,16-17H,1H2,(H,21,22)(H,23,24)(H2,11,12,13)(H2,18,19,20)/p-1/t4-,6-,7-,10-/m1/s1",
            "structure_type": "INCHI",
            "smiles_to_inchi": null,
            "normalized_inchi": "InChI=1S/C10H16N5O13P3/c11-8-5-9(13-2-12-8)15(3-14-5)10-7(17)6(16)4(26-10)1-25-30(21,22)28-31(23,24)27-29(18,19)20/h2-4,6-7,10,16-17H,1H2,(H,21,22)(H,23,24)(H2,11,12,13)(H2,18,19,20)/t4-,6-,7-,10-/m1/s1",
            "canonical_inchi":  "InChI=1S/C10H16N5O13P3/c11-8-5-9(13-2-12-8)15(3-14-5)10-7(17)6(16)4(26-10)1-25-30(21,22)28-31(23,24)27-29(18,19)20/h2-4,6-7,10,16-17H,1H2,(H,21,22)(H,23,24)(H2,11,12,13)(H2,18,19,20)/t4-,6-,7-,10-/m1/s1"
        }
    ]
}
```

For SMILES entries, `smiles_to_inchi` contains the intermediate InChI produced by
conversion. `normalized_inchi` is `null` when no charge layer was removed (only a
SMILES→InChI conversion was performed).

---

## Example files

The `files_examples/` directory contains ready-to-use input files and expected outputs
for all four comparison modes, including InChI files, SMILES files, and MGF files
covering salt normalization, charge normalization, and SMILES-in-MGF cases.

---

## Related repositories

| Repository | Description |
|------------|-------------|
| [inchi-identity-api](https://github.com/alejandraoshea/inchi-identity-api) | Flask REST backend exposing the same comparison engine via HTTP endpoints |
| [inchi-identity-app](https://github.com/alejandraoshea/inchi-identity-app) | Interactive web frontend for browser-based comparison, file upload, and molecular visualization |

---

