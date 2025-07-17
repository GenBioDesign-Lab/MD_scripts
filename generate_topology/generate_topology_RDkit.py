from rdkit import Chem
import os
import argparse

# Parse command line arguments
parser = argparse.ArgumentParser(description='Generate CHARMM RTF topology file from mol2 file')
parser.add_argument('-i', '--input', required=True, help='Input mol2 file')
args = parser.parse_args()

# Get input filename and create output filename
input_file = args.input
output_file = os.path.splitext(input_file)[0] + '.rtf'
pdb_file = os.path.splitext(input_file)[0] + '.pdb'

mol = Chem.MolFromMol2File(input_file, sanitize=False, removeHs=False)

# Check if molecule has more than 1000 atoms for renaming logic
has_many_atoms = mol.GetNumAtoms() > 1000
if has_many_atoms:
    print(f"Molecule has {mol.GetNumAtoms()} atoms. Applying atom renaming")

# Assign custom partial charges based on bonding patterns
def assign_custom_charges(mol):
    """
    Assign custom partial charges based on bonding patterns:
    - C bonded to H: -0.115
    - H bonded to C: 0.115
    - C bonded to O: 0.11
    - O bonded to C: -0.53
    - H bonded to O: 0.42
    - All other C: 0
    """
    for atom in mol.GetAtoms():
        charge = 0.0
        element = atom.GetSymbol()
        
        # Get neighboring atoms
        neighbors = [n.GetSymbol() for n in atom.GetNeighbors()]
        
        if element == 'C':
            if 'H' in neighbors:
                charge = -0.115
            elif 'O' in neighbors:
                charge = 0.11
            else:
                charge = 0.0
        elif element == 'H':
            if 'C' in neighbors:
                charge = 0.115
            elif 'O' in neighbors:
                charge = 0.42
        elif element == 'O':
            if 'C' in neighbors:
                charge = -0.53
        
        atom.SetProp('_CustomCharge', str(charge))

assign_custom_charges(mol)

bonds = [(b.GetBeginAtomIdx(), b.GetEndAtomIdx()) for b in mol.GetBonds()]
angles = []
for b in mol.GetBonds():
    i, j = b.GetBeginAtomIdx(), b.GetEndAtomIdx()
    for n in mol.GetAtomWithIdx(j).GetNeighbors():
        k = n.GetIdx()
        if k != i:
            angles.append((i, j, k))
    for n in mol.GetAtomWithIdx(i).GetNeighbors():
        k = n.GetIdx()
        if k != j:
            angles.append((j, i, k))

angles = list({tuple(a) for a in angles})

torsions = []
for bond12 in mol.GetBonds():
    i, j = bond12.GetBeginAtomIdx(), bond12.GetEndAtomIdx()
    for bond23 in mol.GetAtomWithIdx(j).GetBonds():
        if bond23.GetIdx() == bond12.GetIdx():
            continue
        k = bond23.GetOtherAtomIdx(j)
        for bond34 in mol.GetAtomWithIdx(k).GetBonds():
            if bond34.GetIdx() in (bond12.GetIdx(), bond23.GetIdx()):
                continue
            l = bond34.GetOtherAtomIdx(k)
            torsions.append((i, j, k, l))

torsions = list({t if t < t[::-1] else t[::-1] for t in torsions})

impropers = []
# Try to compute hybridization for improper detection
try:
    mol_copy = Chem.Mol(mol)
    Chem.rdMolOps.SetHybridization(mol_copy)
    
    for i, a in enumerate(mol.GetAtoms()):
        a_copy = mol_copy.GetAtomWithIdx(i)
        nbrs = [n.GetIdx() for n in a.GetNeighbors()]
        if len(nbrs) == 3 and a_copy.GetHybridization() in (Chem.rdchem.HybridizationType.SP2,
                                                           Chem.rdchem.HybridizationType.SP):
            idx, j, k, l = a.GetIdx(), *nbrs
            impropers.append((j, idx, k, l))
except:
    # Fallback: detect impropers based on carbon atoms with 3 aromatic bonds
    for a in mol.GetAtoms():
        if a.GetSymbol() == 'C':
            nbrs = [n.GetIdx() for n in a.GetNeighbors()]
            if len(nbrs) == 3:
                aromatic_bonds = 0
                for bond in a.GetBonds():
                    if bond.GetBondType() == Chem.rdchem.BondType.AROMATIC:
                        aromatic_bonds += 1
                
                if aromatic_bonds >= 2:
                    idx, j, k, l = a.GetIdx(), *nbrs
                    impropers.append((j, idx, k, l))

topology = {
    "bonds": bonds,
    "angles": angles,
    "torsions": torsions,
    "impropers": impropers,
}

# Write CHARMM-style RTF file
def write_charmm_rtf(topology, filename="topology.rtf"):
    """Write topology in CHARMM RTF format"""
    with open(filename, 'w') as f:
        f.write("* Topology file generated from RDKit\n")
        f.write("MASS     1 CCNT     12.011000\n")
        f.write("MASS     2 OCNT     16.000000\n")
        f.write("MASS     3 HOCNT      1.008000\n")
        f.write("MASS     4 HCCNT      1.008000\n")
        f.write("*\n\n")
        # Calculate total charge from custom charges
        total_charge = 0.0
        for i in range(mol.GetNumAtoms()):
            atom = mol.GetAtomWithIdx(i)
            if atom.HasProp('_CustomCharge'):
                charge = float(atom.GetProp('_CustomCharge'))
                total_charge += charge
        
        f.write(f"RESI CNT    {total_charge:8.3f}\n")
        f.write("GROUP\n")
        
        # Atoms section
        atom_counts = {}
        
        for i in range(mol.GetNumAtoms()):
            atom = mol.GetAtomWithIdx(i)
            element = atom.GetSymbol()
            
            if element not in atom_counts:
                atom_counts[element] = 0
            atom_counts[element] += 1
            
            # Apply renaming logic if molecule has >1000 atoms
            count = atom_counts[element]
            if has_many_atoms and count > 999:
                # Use alphabet combinations: C1000->CAA0, C1001->CAA1, etc.
                excess = count - 999
                # Generate all possible 2-letter + 1-digit combinations (AA0-EE9 = 250 combinations)
                letters = 'ABCDE'
                digits = '0123456789'
                combo_index = (excess - 1) % 250  # Cycle through 250 combinations
                
                first_letter_idx = combo_index // 50  # 0-4 for A-E
                second_letter_idx = (combo_index % 50) // 10  # 0-4 for A-E
                digit_idx = combo_index % 10  # 0-9
                
                atom_name = f"{element}{letters[first_letter_idx]}{letters[second_letter_idx]}{digits[digit_idx]}"
            else:
                atom_name = f"{element}{count}"

            # Determine CHARMM type based on actual element (not renamed atom name)
            if element == 'C':
                charmm_type = "CCNT"
            elif element == 'O':
                charmm_type = "OCNT"
            elif element == 'H':
                bonded_atoms = [n for n in atom.GetNeighbors()]
                if bonded_atoms:
                    bonded_element = bonded_atoms[0].GetSymbol()
                    if bonded_element == 'O':
                        charmm_type = "HOCNT"
                    else:
                        charmm_type = "HCCNT"
                else:
                    charmm_type = "HCCNT"
            else:
                charmm_type = element
            
            charge = 0.0
            if atom.HasProp('_CustomCharge'):
                charge = float(atom.GetProp('_CustomCharge'))
            
            f.write(f"ATOM {atom_name:<4s} {charmm_type:<6s} {charge:8.3f}\n")
        f.write("\n")
        
        def get_atom_name(atom_idx):
            atom = mol.GetAtomWithIdx(atom_idx)
            element = atom.GetSymbol()
            count = 1
            for i in range(atom_idx):
                if mol.GetAtomWithIdx(i).GetSymbol() == element:
                    count += 1
            
            # Apply renaming logic if molecule has >1000 atoms
            if has_many_atoms and count > 999:
                # Use alphabet combinations: C1000->CAA0, C1001->CAA1, etc.
                excess = count - 999
                # Generate all possible 2-letter + 1-digit combinations (AA0-EE9 = 250 combinations)
                letters = 'ABCDE'
                digits = '0123456789'
                combo_index = (excess - 1) % 250  # Cycle through 250 combinations
                
                first_letter_idx = combo_index // 50  # 0-4 for A-E
                second_letter_idx = (combo_index % 50) // 10  # 0-4 for A-E
                digit_idx = combo_index % 10  # 0-9
                
                return f"{element}{letters[first_letter_idx]}{letters[second_letter_idx]}{digits[digit_idx]}"
            else:
                return f"{element}{count}"
        
        # Bonds section
        if topology['bonds']:
            for atom1, atom2 in topology['bonds']:
                atom1_name = get_atom_name(atom1)
                atom2_name = get_atom_name(atom2)
                f.write(f"BOND {atom1_name:<4s} {atom2_name:<4s}\n")
            f.write("\n")
        
        # Angles section
        if topology['angles']:
            for atom1, atom2, atom3 in topology['angles']:
                atom1_name = get_atom_name(atom1)
                atom2_name = get_atom_name(atom2)
                atom3_name = get_atom_name(atom3)
                f.write(f"ANGL {atom1_name:<4s} {atom2_name:<4s} {atom3_name:<4s}\n")
            f.write("\n")
        
        # Dihedrals section
        if topology['torsions']:
            for atom1, atom2, atom3, atom4 in topology['torsions']:
                atom1_name = get_atom_name(atom1)
                atom2_name = get_atom_name(atom2)
                atom3_name = get_atom_name(atom3)
                atom4_name = get_atom_name(atom4)
                f.write(f"DIHE {atom1_name:<4s} {atom2_name:<4s} {atom3_name:<4s} {atom4_name:<4s}\n")
            f.write("\n")
        
        # Improper dihedrals section
        if topology['impropers']:
            for atom1, atom2, atom3, atom4 in topology['impropers']:
                atom1_name = get_atom_name(atom1)
                atom2_name = get_atom_name(atom2)
                atom3_name = get_atom_name(atom3)
                atom4_name = get_atom_name(atom4)
                f.write(f"IMPR {atom1_name:<4s} {atom2_name:<4s} {atom3_name:<4s} {atom4_name:<4s}\n")
            f.write("\n")
        
        f.write("END\n")

# Post-process PDB file to rename atoms consistently with RTF topology
def rename_pdb_atoms(pdb_filename, mol):
    """Rename atoms in PDB file to match RTF topology naming (C1, C2, H1, etc.)"""
    def get_atom_name(atom_idx):
        atom = mol.GetAtomWithIdx(atom_idx)
        element = atom.GetSymbol()
        count = 1
        for i in range(atom_idx):
            if mol.GetAtomWithIdx(i).GetSymbol() == element:
                count += 1
        
        # Apply renaming logic if molecule has >1000 atoms
        if has_many_atoms and count > 999:
            # Use alphabet combinations: C1000->CAA0, C1001->CAA1, etc.
            excess = count - 999
            # Generate all possible 2-letter + 1-digit combinations (AA0-EE9 = 250 combinations)
            letters = 'ABCDE'
            digits = '0123456789'
            combo_index = (excess - 1) % 250  # Cycle through 250 combinations
            
            first_letter_idx = combo_index // 50  # 0-4 for A-E
            second_letter_idx = (combo_index % 50) // 10  # 0-4 for A-E
            digit_idx = combo_index % 10  # 0-9
            
            return f"{element}{letters[first_letter_idx]}{letters[second_letter_idx]}{digits[digit_idx]}"
        else:
            return f"{element}{count}"
    
    try:
        with open(pdb_filename, 'r') as f:
            lines = f.readlines()
        
        new_lines = []
        atom_counter = 0
        
        for line in lines:
            if line.startswith('ATOM') or line.startswith('HETATM'):
                if len(line) >= 80:
                    new_atom_name = get_atom_name(atom_counter)
                    new_line = (line[:12] + f"{new_atom_name:<4s}" + line[16:])
                    new_lines.append(new_line)
                    atom_counter += 1
                else:
                    new_lines.append(line)
            else:
                new_lines.append(line)
        
        with open(pdb_filename, 'w') as f:
            f.writelines(new_lines)
        
        return True
    except:
        return False

# Write PDB using OpenBabel
def write_pdb_openbabel(mol2_filename, pdb_filename):
    """Write PDB file using OpenBabel conversion from mol2"""
    try:
        import subprocess
        cmd = ['obabel', mol2_filename, '-O', pdb_filename]
        result = subprocess.run(cmd, capture_output=True, text=True)
        return result.returncode == 0
    except:
        return False

# Alternative: Use Python OpenBabel if available
def write_pdb_openbabel_python(mol2_filename, pdb_filename):
    """Write PDB file using Python OpenBabel library"""
    try:
        from openbabel import openbabel as ob
        conv = ob.OBConversion()
        mol = ob.OBMol()
        conv.SetInAndOutFormats("mol2", "pdb")
        conv.ReadFile(mol, mol2_filename)
        conv.WriteFile(mol, pdb_filename)
        return True
    except:
        return False

# Generate files
write_charmm_rtf(topology, output_file)

pdb_generated = write_pdb_openbabel(input_file, pdb_file)
if not pdb_generated:
    pdb_generated = write_pdb_openbabel_python(input_file, pdb_file)

if pdb_generated:
    rename_pdb_atoms(pdb_file, mol)

print(f"Files generated: {output_file}, {pdb_file}")
print(f"Topology: {len(topology['bonds'])} bonds, {len(topology['angles'])} angles, {len(topology['torsions'])} torsions, {len(topology['impropers'])} impropers")