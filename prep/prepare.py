"""Extract 4EY7 chain A and prepare the E20 ligand and Vina configuration."""

from pathlib import Path

import gemmi
import numpy as np
from meeko import MoleculePreparation, PDBQTWriterLegacy
from rdkit import Chem


ROOT = Path(__file__).resolve().parents[1]
DOCK = ROOT / "docking"
DATA = ROOT / "data"


def xyz(line):
    return [float(line[start : start + 8]) for start in (30, 38, 46)]


def write_pdb(path, lines):
    path.write_text("\n".join([*lines, "TER", "END"]) + "\n", encoding="ascii")


def extract_chain():
    rec, ligand = [], {}
    original = (DATA / "4EY7.pdb").read_text(encoding="ascii")
    disulfides = [line for line in original.splitlines()
                  if line.startswith("SSBOND") and line[15] == line[29] == "A"]
    for line in original.splitlines():
        if line[:6] not in {"ATOM  ", "HETATM"} or line[21] != "A":
            continue
        if line[16] not in {" ", "A"}:
            continue
        line = line[:16] + " " + line[17:]
        if line.startswith("ATOM  "):
            rec.append(line)
        elif line[17:26] == "E20 A 604" and line[76:78].strip() != "H":
            name = line[12:16].strip()
            if name in ligand:
                raise ValueError(f"Repeated crystallographic atom {name}")
            ligand[name] = line
    if len(rec) != 4130 or len(ligand) != 28 or len(disulfides) != 3:
        raise ValueError("Unexpected atom or disulfide count in the 4EY7 input")
    write_pdb(DOCK / "receptor_A.pdb", [*disulfides, *rec])
    return ligand


def reference_molecule(ligand):
    ccd = gemmi.cif.read_file(str(DATA / "E20.cif")).sole_block()
    graph, indexes = Chem.RWMol(), {}
    for row in ccd.find("_chem_comp_atom.", ["atom_id", "type_symbol", "charge"]):
        name, element, charge = list(row)
        if element == "H":
            continue
        atom = Chem.Atom(element)
        atom.SetFormalCharge(int(charge))
        indexes[name] = graph.AddAtom(atom)
    if set(indexes) != set(ligand):
        raise ValueError("E20 CCD and 4EY7 chain A atom names disagree")
    bond_type = {"SING": Chem.BondType.SINGLE, "DOUB": Chem.BondType.DOUBLE}
    for row in ccd.find("_chem_comp_bond.", ["atom_id_1", "atom_id_2", "value_order"]):
        first, second, order = list(row)
        if first in indexes and second in indexes:
            graph.AddBond(indexes[first], indexes[second], bond_type[order])
    conformer = Chem.Conformer(len(indexes))
    for name, position in indexes.items():
        conformer.SetAtomPosition(position, xyz(ligand[name]))
    graph.AddConformer(conformer)
    molecule = graph.GetMol()
    Chem.SanitizeMol(molecule)
    Chem.AssignStereochemistryFrom3D(molecule)
    return molecule, indexes


def write_sdf(path, molecule):
    with Chem.SDWriter(str(path)) as writer:
        writer.write(molecule)


def main():
    ligand = extract_chain()
    reference, indexes = reference_molecule(ligand)
    write_sdf(DOCK / "E20_crystal.sdf", reference)
    # Keep original PDB atom names for the contact figure, with CCD connectivity.
    serials = {index: int(ligand[name][6:11]) for name, index in indexes.items()}
    connectivity = [
        f"CONECT{serials[bond.GetBeginAtomIdx()]:5d}{serials[bond.GetEndAtomIdx()]:5d}"
        for bond in reference.GetBonds()
    ]
    write_pdb(DOCK / "E20_crystal.pdb", [*ligand.values(), *connectivity])

    # The docking model uses a protonated tertiary piperidine nitrogen.
    docking_mol = Chem.Mol(reference)
    docking_mol.GetAtomWithIdx(indexes["N14"]).SetFormalCharge(1)
    Chem.SanitizeMol(docking_mol)
    docking_mol = Chem.AddHs(docking_mol, addCoords=True)
    write_sdf(DOCK / "E20_docking.sdf", docking_mol)
    (setup,) = MoleculePreparation().prepare(docking_mol)
    pdbqt, valid, message = PDBQTWriterLegacy.write_string(setup)
    if not valid:
        raise ValueError(message)
    (DOCK / "E20_docking.pdbqt").write_text(pdbqt, encoding="ascii")

    coordinates = reference.GetConformer().GetPositions()
    center = coordinates.mean(axis=0)
    diameter = np.linalg.norm(coordinates[:, None, :] - coordinates[None, :, :], axis=2).max()
    config = [
        "receptor = docking/receptor_A.pdbqt",
        "ligand = docking/E20_docking.pdbqt",
        *(f"center_{axis} = {value:.4f}" for axis, value in zip("xyz", center)),
        *(f"size_{axis} = 24" for axis in "xyz"),
        "exhaustiveness = 32",
        "num_modes = 9",
        "energy_range = 4",
        "seed = 2027",
        "cpu = 4",
        "out = docking/vina_poses.pdbqt",
    ]
    (DOCK / "vina_config.txt").write_text("\n".join(config) + "\n", encoding="ascii")
    print(f"E20 center {center.round(4)}, heavy-atom diameter {diameter:.2f} angstrom")
    print(f"Crystal C8 configuration: {reference.GetAtomWithIdx(indexes['C8']).GetProp('_CIPCode')}")


if __name__ == "__main__":
    main()
