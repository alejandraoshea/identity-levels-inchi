# inchi-identity

A hierarchical molecular identity comparison framework for untargeted metabolomics,
implemented as a Python CLI and published on PyPI.

Instead of binary exact matching, this tool evaluates molecular equivalence across six
progressive normalization layers, returning a structured equivalence profile that reflects
the structural resolution actually supported by the experimental evidence. For lipids, a
dedicated four-level hierarchy (Levels A–D) addresses cis/trans geometry, sn-position,
intra-chain double bond position, and global sum composition — each reported individually.

Both InChI strings and SMILES are accepted as input — SMILES are automatically converted
to canonical InChI before comparison.

---

## Comparison Layers

| Layer | Key | Name | Description |
|-------|-----|------|-------------|
| 1 | `complete_identity` | Complete Identity | Exact InChI string equality |
| 2 | `isotope` | Isotopic Independence | Equality after `/i` layer removal |
| 3 | `salt` | Salt Independence | Equality after counterion removal (RDKit SaltRemover) |
| 4 | `charge` | Charge Independence | Equality after charge/protonation normalization |
| 5 | `double_bond` | Double Bond Independence | Position-independent comparison (any of A–D match) |
| 5A | `isomer_level_a` | (A) Cis/Trans | Removes `/b` cis/trans geometry; applies to all molecules |
| 5B | `isomer_level_b` | (B) sn-Position | Compares acyl chains as an unordered set; lipids only |
| 5C | `isomer_level_c` | (C) Intra-Chain Position | Discards double bond position within each chain; lipids only |
| 5D | `isomer_level_d` | (D) Sum Composition | Total carbon and double bond count only; lipids only |
| 6 | `tautomer` | Tautomeric Independence | Equality after canonical tautomer generation |

Preprocessing is **cumulative**: each layer applies all normalizations from preceding
layers before performing its own comparison.

### Lipid sub-levels (Layer 5A–D)

Layer 5 (`double_bond`) returns a single `true/false` result: `true` if the two molecules
match at *any* of the four sub-levels. Layers 5A–D are reported **individually** — each one
tells you at which level of structural resolution the pair first becomes equivalent.

Sub-levels are also **cumulative downward**: if a pair matches at Level A, it
automatically matches at B, C, and D as well (since each is more lenient).

| Level | Key | Match criterion |
|-------|-----|----------------|
| A | `isomer_level_a` | Identical after removing cis/trans geometry |
| B | `isomer_level_b` | Same chains regardless of which sn-position they occupy |
| C | `isomer_level_c` | Same per-chain composition (C, DB, O) ignoring double bond position |
| D | `isomer_level_d` | Same total carbon count and total double bond count |

For non-lipids, levels B, C, and D are not applicable — they are absent from the result
(`N/A` in the web interface).

Lipid classification uses the ClassyFire API (5 s timeout, InChIKey lookup), with a
SMARTS headgroup pattern fallback when ClassyFire is unavailable or the compound is not
in its database. Headgroup validation uses a library of 428 SMARTS patterns (79 primary +
255 generated at runtime by combining 15 glycolipid templates with 17 monosaccharide
patterns + 94 sugar-agnostic backbone patterns).

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

Runs every layer enabled in `default_config.json` (all layers on by default):

```bash
inchi compare-pair "<inchi_1>" "<inchi_2>"
inchi compare-pair "<inchi_1>" "<inchi_2>" 2>/dev/null > result.json
```

SMILES are also accepted:

```bash
inchi compare-pair "CC(=O)O" "CC(=O)[O-].[Na+]"
```

#### Using a custom config file

To run only a subset of layers without editing the default config, create a custom JSON
file and pass it with `--config`:

```bash
inchi compare-pair "<inchi_1>" "<inchi_2>" --config path/to/my_config.json
```

A config file has the same structure as `default_config.json` — set any layer to `false`
to skip it. Example that runs only the four isomer sub-levels:

```json
{
  "identity_criteria": {
    "complete_identity":   { "enabled": false },
    "isotope_independence":  { "isotope_independent_identity": false },
    "salt_independence":     { "desalted_identity": false },
    "charge_independence":   { "charge_independent_identity": false },
    "isomer_independence": {
      "double_bond_position_independent_identity": false,
      "isomer_level_a_cis_trans":      true,
      "isomer_level_b_sn_position":    true,
      "isomer_level_c_intra_chain":    true,
      "isomer_level_d_sum_composition": true
    },
    "tautomer_independence": { "tautomer_independent_identity": false, "inchitrust_path": null }
  }
}
```

A ready-to-use version of this file is included at
`src/inchi_identity/configs/isomer_only_config.json`.

Redirect output to a file (use `2>/dev/null` on Linux/macOS or `2>$null` on PowerShell
to suppress the loading messages from stderr):

```bash
# Linux / macOS
inchi compare-pair "<inchi_1>" "<inchi_2>" \
    --config src/inchi_identity/configs/isomer_only_config.json \
    2>/dev/null > result.json

# PowerShell
inchi compare-pair "<inchi_1>" "<inchi_2>" `
    --config src/inchi_identity/configs/isomer_only_config.json `
    2>$null > result.json
```

### Mode 2 — Compare with selected layers only (no config editing needed)

Pass layer keys directly on the command line with `--layers`. Only the specified layers
are evaluated; the config file is ignored for layer selection.

```bash
inchi compare-pair-layers "<inchi_1>" "<inchi_2>" \
    --layers isotope charge double_bond tautomer
```

**All available layer keys:**

| Key | Description |
|-----|-------------|
| `complete_identity` | Exact match |
| `isotope` | Isotope independence |
| `salt` | Salt independence |
| `charge` | Charge independence |
| `double_bond` | Full double bond / stereochemical independence (A–D combined) |
| `isomer_level_a` | (A) Cis/Trans |
| `isomer_level_b` | (B) sn-Position |
| `isomer_level_c` | (C) Intra-Chain Position |
| `isomer_level_d` | (D) Sum Composition |
| `tautomer` | Tautomer independence |

Multiple keys can be combined:

```bash
# Only the four isomer sub-levels:
inchi compare-pair-layers "<inchi_1>" "<inchi_2>" \
    --layers isomer_level_a isomer_level_b isomer_level_c isomer_level_d

# Double bond + all four sub-levels reported separately:
inchi compare-pair-layers "<inchi_1>" "<inchi_2>" \
    --layers double_bond isomer_level_a isomer_level_b isomer_level_c isomer_level_d

# Classic layers only:
inchi compare-pair-layers "<inchi_1>" "<inchi_2>" \
    --layers complete_identity isotope salt charge tautomer
```

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

A custom config can also be passed to file-based comparison:

```bash
inchi compare file1.txt file2.txt --config my_config.json
```

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
ISOMER_INDEPENDENCE_A
ISOMER_INDEPENDENCE_B
ISOMER_INDEPENDENCE_C
ISOMER_INDEPENDENCE_D
TAUTOMER_INDEPENDENCE
```

The pipeline: (1) parses both files and extracts identifiers, (2) internally normalizes
each file so equivalent entries share a canonical InChI, (3) performs a cross-file
comparison and replaces identifiers in file 2 with the canonical form from file 1,
(4) writes the unified MGF and a JSON log recording every transformation. Entries without
identifiers and all spectral peak data are preserved unchanged.

---

## Output format

### compare-pair — all layers

```json
{
    "inchi_1": "InChI=1S/C18H34O2/.../b10-9-",
    "inchi_2": "InChI=1S/C18H34O2/.../b10-9+",
    "results": {
        "COMPLETE_IDENTITY": false,
        "ISOTOPIC_INDEPENDENCE": false,
        "SALTS_INDEPENDENCE": false,
        "CHARGES_INDEPENDENCE": false,
        "DOUBLE_BONDS_INDEPENDENCE": true,
        "ISOMER_INDEPENDENCE_A": true,
        "ISOMER_INDEPENDENCE_B": true,
        "ISOMER_INDEPENDENCE_C": true,
        "ISOMER_INDEPENDENCE_D": true,
        "TAUTOMER_INDEPENDENCE": true
    }
}
```

### compare-pair — isomer sub-levels only (different double bond position)

cis-Oleic acid (Δ9) vs Vaccenic acid (Δ11): same 18:1 composition, different position.

```json
{
    "inchi_1": "InChI=1S/C18H34O2/.../b10-9-",
    "inchi_2": "InChI=1S/C18H34O2/.../b12-11+",
    "results": {
        "ISOMER_INDEPENDENCE_A": false,
        "ISOMER_INDEPENDENCE_B": false,
        "ISOMER_INDEPENDENCE_C": true,
        "ISOMER_INDEPENDENCE_D": true
    }
}
```

Level A and B fail (position differs); Level C matches (both are 18:1 per chain);
Level D matches (total C=18, DB=1).

### compare-pair — charge normalization

```json
{
    "inchi_1": "InChI=1S/C5H11NO2/c1-6(2,3)4-5(7)8/h4H2,1-3H3/p+1",
    "inchi_2": "InChI=1S/C5H11NO2/c1-6(2,3)4-5(7)8/h4H2,1-3H3",
    "results": {
        "COMPLETE_IDENTITY": false,
        "ISOTOPIC_INDEPENDENCE": false,
        "SALTS_INDEPENDENCE": false,
        "CHARGES_INDEPENDENCE": true,
        "DOUBLE_BONDS_INDEPENDENCE": true,
        "ISOMER_INDEPENDENCE_A": true,
        "ISOMER_INDEPENDENCE_B": null,
        "ISOMER_INDEPENDENCE_C": null,
        "ISOMER_INDEPENDENCE_D": null,
        "TAUTOMER_INDEPENDENCE": true
    }
}
```

`null` indicates the level is not applicable (non-lipid molecule).

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
            "original_structure": "InChI=1S/.../p-1/t4-,6-,7-,10-/m1/s1",
            "structure_type": "INCHI",
            "smiles_to_inchi": null,
            "normalized_inchi": "InChI=1S/.../t4-,6-,7-,10-/m1/s1",
            "canonical_inchi":  "InChI=1S/.../t4-,6-,7-,10-/m1/s1"
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
