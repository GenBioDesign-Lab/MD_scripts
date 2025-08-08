import argparse
import os
import subprocess
from typing import Dict, List, Tuple

from rdkit import Chem


# ------------------------------
# CLI
# ------------------------------
def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Generate CHARMM RTF topology from a MOL2 file and a PDB via OpenBabel"
    )
    parser.add_argument("-i", "--input", required=True, help="Input .mol2 file")
    parser.add_argument(
        "-r", "--resname", default="CNT", help="Residue name to write in RTF and PDB"
    )
    parser.add_argument(
        "--pdb-only", action="store_true", help="Only generate PDB (skip RTF)"
    )
    parser.add_argument(
        "--rtf-only", action="store_true", help="Only generate RTF (skip PDB)"
    )
    return parser.parse_args()


# ------------------------------
# Atom naming
# ------------------------------
def base36(n: int) -> str:
    digits = "0123456789ABCDEFGHIJKLMNOPQRSTUVWXYZ"
    if n == 0:
        return "0"
    out = []
    while n > 0:
        n, r = divmod(n, 36)
        out.append(digits[r])
    return "".join(reversed(out))


def build_atom_name_map(mol: Chem.Mol) -> List[str]:
    """
    Produce a stable 4-character atom name for each atom index, consistent across RTF and PDB.

    Strategy:
    - Per element counter in traversal order
    - Name = <ElementSymbol><3-char base36 counter>, padded to width 3
      e.g. C -> C001.. C00A .. C010, H -> H001..
    - Ensures max width 4 for PDB atom name column
    """
    element_counts: Dict[str, int] = {}
    names: List[str] = []
    for atom in mol.GetAtoms():
        element = atom.GetSymbol()
        count = element_counts.get(element, 0) + 1
        element_counts[element] = count
        suffix = base36(count).rjust(3, "0")  # width 3, base36
        name = f"{element}{suffix}"[:4]  # hard cap to 4 chars for PDB
        names.append(name)
    return names


# ------------------------------
# Custom charges
# ------------------------------
def assign_custom_charges(mol: Chem.Mol) -> None:
    for atom in mol.GetAtoms():
        element = atom.GetSymbol()
        neighbor_symbols = [n.GetSymbol() for n in atom.GetNeighbors()]

        charge = 0.0
        if element == "C":
            if "H" in neighbor_symbols:
                charge = -0.115
            elif "O" in neighbor_symbols:
                charge = 0.11
        elif element == "H":
            if "C" in neighbor_symbols:
                charge = 0.115
            elif "O" in neighbor_symbols:
                charge = 0.42
        elif element == "O":
            if "C" in neighbor_symbols:
                charge = -0.53

        atom.SetProp("_CustomCharge", str(charge))


# ------------------------------
# Topology extraction (bonds / angles / torsions / impropers)
# ------------------------------
def extract_bonds(mol: Chem.Mol) -> List[Tuple[int, int]]:
    return [(b.GetBeginAtomIdx(), b.GetEndAtomIdx()) for b in mol.GetBonds()]


def extract_angles(mol: Chem.Mol) -> List[Tuple[int, int, int]]:
    angles: List[Tuple[int, int, int]] = []
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
    return list({tuple(a) for a in angles})


def extract_torsions(mol: Chem.Mol) -> List[Tuple[int, int, int, int]]:
    torsions: List[Tuple[int, int, int, int]] = []
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
    return list({t if t < t[::-1] else t[::-1] for t in torsions})


def extract_impropers(mol: Chem.Mol) -> List[Tuple[int, int, int, int]]:
    impropers: List[Tuple[int, int, int, int]] = []
    try:
        mol_copy = Chem.Mol(mol)
        Chem.rdMolOps.SetHybridization(mol_copy)
        for i, a in enumerate(mol.GetAtoms()):
            a_copy = mol_copy.GetAtomWithIdx(i)
            nbrs = [n.GetIdx() for n in a.GetNeighbors()]
            if len(nbrs) == 3 and a_copy.GetHybridization() in (
                Chem.rdchem.HybridizationType.SP2,
                Chem.rdchem.HybridizationType.SP,
            ):
                idx, j, k, l = a.GetIdx(), *nbrs
                impropers.append((j, idx, k, l))
    except Exception:
        for a in mol.GetAtoms():
            if a.GetSymbol() == "C":
                nbrs = [n.GetIdx() for n in a.GetNeighbors()]
                if len(nbrs) == 3:
                    aromatic_bonds = sum(
                        1 for b in a.GetBonds() if b.GetBondType() == Chem.rdchem.BondType.AROMATIC
                    )
                    if aromatic_bonds >= 2:
                        idx, j, k, l = a.GetIdx(), *nbrs
                        impropers.append((j, idx, k, l))
    return impropers


# ------------------------------
# RTF writer
# ------------------------------
def write_charmm_rtf(
    mol: Chem.Mol,
    atom_names: List[str],
    bonds: List[Tuple[int, int]],
    angles: List[Tuple[int, int, int]],
    torsions: List[Tuple[int, int, int, int]],
    impropers: List[Tuple[int, int, int, int]],
    filename: str,
    residue_name: str,
) -> None:
    total_charge = 0.0
    for i in range(mol.GetNumAtoms()):
        atom = mol.GetAtomWithIdx(i)
        if atom.HasProp("_CustomCharge"):
            total_charge += float(atom.GetProp("_CustomCharge"))

    with open(filename, "w") as f:
        f.write("* Topology file generated from RDKit\n")
        f.write("MASS     1 CCNT     12.011000\n")
        f.write("MASS     2 OCNT     16.000000\n")
        f.write("MASS     3 HOCNT      1.008000\n")
        f.write("MASS     4 HCCNT      1.008000\n")
        f.write("*\n\n")
        f.write(f"RESI {residue_name}    {total_charge:8.3f}\n")
        f.write("GROUP\n")

        # atoms
        for i in range(mol.GetNumAtoms()):
            atom = mol.GetAtomWithIdx(i)
            element = atom.GetSymbol()
            if element == "C":
                charmm_type = "CCNT"
            elif element == "O":
                charmm_type = "OCNT"
            elif element == "H":
                bonded_atoms = [n for n in atom.GetNeighbors()]
                if bonded_atoms and bonded_atoms[0].GetSymbol() == "O":
                    charmm_type = "HOCNT"
                else:
                    charmm_type = "HCCNT"
            else:
                charmm_type = element

            charge = float(atom.GetProp("_CustomCharge")) if atom.HasProp("_CustomCharge") else 0.0
            f.write(f"ATOM {atom_names[i]:<4s} {charmm_type:<6s} {charge:8.3f}\n")
        f.write("\n")

        def nm(idx: int) -> str:
            return atom_names[idx]

        # bonds
        for i, j in bonds:
            f.write(f"BOND {nm(i):<4s} {nm(j):<4s}\n")
        f.write("\n")

        # angles
        for i, j, k in angles:
            f.write(f"ANGL {nm(i):<4s} {nm(j):<4s} {nm(k):<4s}\n")
        f.write("\n")

        # torsions
        for i, j, k, l in torsions:
            f.write(f"DIHE {nm(i):<4s} {nm(j):<4s} {nm(k):<4s} {nm(l):<4s}\n")
        f.write("\n")

        # impropers
        for i, j, k, l in impropers:
            f.write(f"IMPR {nm(i):<4s} {nm(j):<4s} {nm(k):<4s} {nm(l):<4s}\n")
        f.write("\nEND\n")


# ------------------------------
# PDB via OpenBabel CLI, then rename atom names
# ------------------------------
def write_pdb_with_openbabel_cli(mol2_filename: str, pdb_filename: str) -> None:
    cmd = ["obabel", mol2_filename, "-O", pdb_filename]
    result = subprocess.run(cmd, capture_output=True, text=True)
    if result.returncode != 0:
        raise RuntimeError(f"OpenBabel failed: {result.stderr.strip()}")


def rename_pdb_atoms_inplace(pdb_filename: str, atom_names: List[str], residue_name: str) -> None:
    with open(pdb_filename, "r") as f:
        lines = f.readlines()

    new_lines: List[str] = []
    atom_counter = 0
    for line in lines:
        if line.startswith("ATOM") or line.startswith("HETATM"):
            if len(line) < 54:
                new_lines.append(line)
                continue
            # PDB columns: atom name (13-16), res name (18-20)
            atom_name = f"{atom_names[atom_counter]:<4s}"[:4]
            line = line[:12] + atom_name + line[16:]
            # set residue name to requested value
            if len(line) >= 20:
                line = line[:17] + f"{residue_name:>3s}" + line[20:]
            new_lines.append(line)
            atom_counter += 1
        else:
            new_lines.append(line)

    with open(pdb_filename, "w") as f:
        f.writelines(new_lines)


# ------------------------------
# Main
# ------------------------------
def main() -> None:
    args = parse_args()
    input_file = os.path.abspath(args.input)
    rtf_file = os.path.splitext(input_file)[0] + ".rtf"
    pdb_file = os.path.splitext(input_file)[0] + ".pdb"

    mol = Chem.MolFromMol2File(input_file, sanitize=False, removeHs=False)
    if mol is None:
        raise RuntimeError(f"Failed to read MOL2: {input_file}")

    # naming map up-front for consistent usage
    atom_names = build_atom_name_map(mol)

    # custom charges
    assign_custom_charges(mol)

    # topology
    bonds = extract_bonds(mol)
    angles = extract_angles(mol)
    torsions = extract_torsions(mol)
    impropers = extract_impropers(mol)

    # outputs
    if not args.pdb_only:
        write_charmm_rtf(
            mol,
            atom_names,
            bonds,
            angles,
            torsions,
            impropers,
            rtf_file,
            residue_name=args.resname,
        )

    if not args.rtf_only:
        write_pdb_with_openbabel_cli(input_file, pdb_file)
        rename_pdb_atoms_inplace(pdb_file, atom_names, residue_name=args.resname)

    print(
        f"Files generated: {rtf_file if not args.pdb_only else '(skipped)'}, "
        f"{pdb_file if not args.rtf_only else '(skipped)'}"
    )
    print(
        "Topology: "
        f"{len(bonds)} bonds, {len(angles)} angles, {len(torsions)} torsions, {len(impropers)} impropers"
    )


if __name__ == "__main__":
    main()


