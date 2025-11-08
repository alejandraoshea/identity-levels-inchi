from rdkit import Chem
from rdkit.Chem import Descriptors
from rdkit.Chem import Draw

smiles = 'CC(=O)OC1=CC=CC=C1C(=O)O'
molecula = Chem.MolFromSmiles(smiles)
Draw.MolToImage(molecula).show()

peso_molecular = Descriptors.MolWt(molecula)
print(f"Masa: {peso_molecular:.3f} g/mol")
# Masa: 180.159 g/mol

methane = Chem.MolFromSmiles("C")
methane

# this phenylalanine uses canonical smiles (no stereochemistry information)
phenylalanine = Chem.MolFromSmiles("C1=CC=C(C=C1)CC(C(=O)O)N")
phenylalanine

# phenyalanine using isomeric smiles
l_phenylalanine = Chem.MolFromSmiles("C1=CC=C(C=C1)C[C@@H](C(=O)O)N")
l_phenylalanine

benzene = Chem.MolFromSmiles("c1ccccc1")
benzene

num_methane = methane.GetNumAtoms()
print(f"The number of atoms in methane is {num_methane}")

num_methane_h = methane.GetNumAtoms(onlyExplicit=False)
print(f"The number of atoms in methane including hydrogens is {num_methane_h}.")

ethanol = Chem.MolFromSmiles("CCO")

for atom in ethanol.GetAtoms():
    print(atom.GetSymbol(), atom.GetMass())

ethanol_h = Chem.AddHs(ethanol)
ethanol_h

for atom in ethanol_h.GetAtoms():
    print(atom.GetSymbol())


for bond in ethanol.GetBonds():
    print(bond.GetBondType())

print("Printing info for benzene:")
print(f"The molecular weight is {Descriptors.MolWt(benzene)}")
print(f"The number of aromatic rings is {Descriptors.NumAromaticRings(benzene)}")