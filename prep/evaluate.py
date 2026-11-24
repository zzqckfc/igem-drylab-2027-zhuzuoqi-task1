"""Measure heavy-atom RMSD in the receptor coordinate frame."""

import csv
import json
import re

import numpy as np
from rdkit import Chem

from prepare import DOCK


def measured_rmsd(reference, pose):
    # Protonation is a docking setting. Normalize N charge for graph comparison.
    pose = Chem.Mol(pose)
    for atom in pose.GetAtoms():
        if atom.GetAtomicNum() == 7 and atom.GetFormalCharge() == 1:
            atom.SetFormalCharge(0)
            atom.SetNumExplicitHs(0)
    Chem.SanitizeMol(pose)
    if Chem.MolToSmiles(reference) != Chem.MolToSmiles(pose):
        raise ValueError("Pose and reference have differing chemistry or stereochemistry")
    mappings = pose.GetSubstructMatches(reference, uniquify=False, useChirality=True)
    if not mappings:
        raise ValueError("Unable to map E20 heavy atoms")
    crystal = reference.GetConformer().GetPositions()
    predicted = pose.GetConformer().GetPositions()
    distances = [
        np.sqrt(np.mean(np.sum((crystal - predicted[list(mapping)]) ** 2, axis=1)))
        for mapping in mappings
    ]
    return float(min(distances)), len(mappings)


def main():
    (reference,) = Chem.SDMolSupplier(str(DOCK / "E20_crystal.sdf"))
    poses = list(Chem.SDMolSupplier(str(DOCK / "vina_poses.sdf")))
    scores = re.findall(r"^REMARK VINA RESULT:\s+(-?\d+\.\d+)",
                        (DOCK / "vina_poses.pdbqt").read_text(), re.MULTILINE)
    if not poses or any(pose is None for pose in poses) or len(poses) != len(scores):
        raise ValueError("Mismatch between Vina scores and exported SDF poses")

    with (DOCK / "pose_metrics.tsv").open("w", encoding="ascii", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t")
        writer.writerow(["rank", "vina_score_kcal_mol", "crystal_heavy_atom_rmsd_angstrom"])
        for rank, (pose, score) in enumerate(zip(poses, scores), 1):
            rmsd, mappings = measured_rmsd(reference, pose)
            writer.writerow([rank, score, f"{rmsd:.4f}"])
            if rank == 1:
                summary = {
                    "pdb_id": "4EY7",
                    "reference": "A:E20:604",
                    "rank": rank,
                    "vina_score_kcal_mol": float(score),
                    "crystal_heavy_atom_rmsd_angstrom": round(rmsd, 4),
                    "heavy_atoms": reference.GetNumAtoms(),
                    "symmetry_mappings": mappings,
                    "success_threshold_angstrom": 2,
                    "success": bool(rmsd < 2.0),
                }
    (DOCK / "result.json").write_text(json.dumps(summary, indent=2) + "\n", encoding="ascii")
    with Chem.SDWriter(str(DOCK / "vina_rank1.sdf")) as writer:
        writer.write(poses[0])
    print(json.dumps(summary, indent=2))


if __name__ == "__main__":
    main()
