from rdkit import Chem
from rdkit.Chem import Draw
#draw the structure of the molecule from mol2 file
mol = Chem.MolFromMol2File("/data01/genbiolab/mdanh/data/MD_scripts/test/cnt6/cnt6_45/raw/cnt_6_45_oh.mol2", sanitize=False)
Draw.MolToFile(mol, "structure.png", size=(3000,3000))