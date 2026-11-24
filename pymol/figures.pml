# Run from the repository root: python -m pymol -cq pymol/figures.pml
python
import json
from pathlib import Path

import numpy as np
from pymol import cmd

OUT = Path("pymol/_render")
OUT.mkdir(parents=True, exist_ok=True)

cmd.reinitialize()
cmd.load("docking/receptor_A.pdb", "ache")
cmd.load("docking/E20_crystal.pdb", "crystal")
cmd.load("docking/vina_rank1.sdf", "rank1")
cmd.remove("hydro")
cmd.hide("everything", "all")
cmd.bg_color("white")
cmd.set("ray_opaque_background", 0)
cmd.set("orthoscopic", 1)
cmd.set("antialias", 2)
cmd.set("ray_shadow", 0)
cmd.set("two_sided_lighting", 1)
cmd.set("ambient", 0.50)
cmd.set("direct", 0.58)
cmd.set("specular", 0.24)
cmd.set("shininess", 26)
cmd.set("ray_trace_fog", 0)
cmd.set("depth_cue", 0)
cmd.set("valence", 0)
cmd.set("cartoon_fancy_helices", 1)
cmd.set("max_threads", 6)

colors = {
    "enzyme": [0.60, 0.69, 0.77],
    "gorge": [0.33, 0.51, 0.62],
    "surface_gray": [0.59, 0.71, 0.79],
    "residue": [0.38, 0.49, 0.57],
    "experimental": [0.00, 0.48, 0.46],
    "predicted": [0.88, 0.36, 0.11],
    "oxygen": [0.78, 0.17, 0.29],
    "nitrogen": [0.20, 0.35, 0.77],
    "hydrogen_bond": [0.59, 0.30, 0.62],
}
for name, color in colors.items():
    cmd.set_color(name, color)

crystal_coords = cmd.get_coords("crystal")
center = crystal_coords.mean(axis=0)


def unit(vector):
    return vector / np.linalg.norm(vector)


def camera(focus, long_axis, vertical_hint, half_height, view_distance, depth, roll=0,
           turn=0):
    x_axis = unit(long_axis)
    y_axis = unit(vertical_hint - x_axis * np.dot(vertical_hint, x_axis))
    z_axis = np.cross(x_axis, y_axis)
    turn_radians = np.deg2rad(turn)
    y_axis = np.cos(turn_radians)*y_axis + np.sin(turn_radians)*z_axis
    z_axis = np.cross(x_axis, y_axis)
    angle = np.deg2rad(roll)
    cos, sin = np.cos(angle), np.sin(angle)
    rotation = np.array([cos*x_axis - sin*y_axis, sin*x_axis + cos*y_axis, z_axis])
    fov = np.rad2deg(2*np.arctan(half_height/view_distance))
    cmd.set("field_of_view", float(fov))
    cmd.set_view([*rotation.T.ravel(), 0, 0, -view_distance, *focus,
                  view_distance-depth, view_distance+depth, 1])
    return rotation


def ray(filename, width, height):
    cmd.png(str(OUT / filename), width, height, dpi=300, ray=1, quiet=1)


def point(selection):
    return np.asarray(cmd.get_atom_coords(selection))


def project(selection, focus, rotation, half_height, width, height):
    xyz = np.mean(cmd.get_coords(selection), axis=0) - focus
    x, y, _ = rotation @ xyz
    scale = height / (2 * half_height)
    return [float(width/2 + scale*x), float(height/2 - scale*y)]


def ligand_element_colors():
    cmd.color("experimental", "crystal")
    cmd.color("oxygen", "crystal and elem O")
    cmd.color("nitrogen", "crystal and elem N")


# Figure 1. The same orientation connects the full fold to the pocket surface.
cmd.color("enzyme", "ache")
cmd.color("gorge", "byres (ache within 6 of crystal)")
cmd.show("cartoon", "ache")
cmd.show("sticks", "crystal")
cmd.show("spheres", "crystal")
cmd.set("sphere_scale", 0.33, "crystal")
cmd.set("stick_radius", 0.33, "crystal")
ligand_element_colors()
overview_center = np.mean(cmd.get_extent("ache"), axis=0)
overview_x = np.array([0.88, -0.24, 0.42])
opening = center - cmd.get_coords("ache and name CA").mean(axis=0)
overview_y = np.cross(opening, overview_x)
full_view = camera(overview_center, overview_x, overview_y, 33, 130, 60, roll=-8)
ray("overview.png", 1210, 1130)
overview_anchor = project("crystal and name C26", overview_center, full_view, 33, 1210, 1130)

cmd.hide("cartoon", "ache")
cmd.hide("spheres", "crystal")
cmd.set("surface_carve_selection", "crystal", "ache")
cmd.set("surface_carve_cutoff", 6.0, "ache")
cmd.set("surface_color", "surface_gray", "ache")
cmd.set("transparency", 0.55, "ache")
cmd.set("surface_quality", 2, "ache")
cmd.show("surface", "ache")
cmd.set("stick_radius", 0.29, "crystal")
camera(center, overview_x, overview_y, 13.0, 62, 26, roll=-8)
ray("gorge.png", 1210, 1130)

# Figure 2. Compare two projections at the same scale, 90 degrees apart.
long_axis = point("crystal and name C21") - point("crystal and name C2")
vertical = point("ache and resi 295 and name N") - center
cmd.hide("everything", "all")
cmd.show("sticks", "crystal or rank1")
cmd.color("experimental", "crystal")
cmd.color("predicted", "rank1")
cmd.set("stick_radius", 0.20, "crystal")
cmd.set("stick_radius", 0.16, "rank1")
cmd.set("stick_transparency", 0.12, "crystal")
camera(center, long_axis, vertical, 5.35, 38, 23, roll=10)
ray("overlay_a.png", 1210, 720)
camera(center, long_axis, vertical, 5.35, 38, 23, roll=10, turn=90)
ray("overlay_b.png", 1210, 720)

# Figure 3. Selected residue side chains and the backbone donor of Phe295.
cmd.hide("everything", "all")
cmd.select("key", "ache and resi 86+286+295+337")
cmd.color("residue", "key")
cmd.show("sticks", "(key and not name C+O+N) or (key and resi 295 and name N)")
cmd.set("stick_radius", 0.17, "key")
cmd.color("oxygen", "key and elem O")
cmd.color("nitrogen", "key and elem N")
cmd.show("sticks", "crystal")
cmd.set("stick_radius", 0.27, "crystal")
cmd.set("stick_transparency", 0, "crystal")
ligand_element_colors()
hb_distance = cmd.distance("hbond", "crystal and name O24",
                            "ache and resi 295 and name N")
cmd.hide("labels", "hbond")
cmd.set("dash_color", "hydrogen_bond", "hbond")
cmd.set("dash_radius", 0.070, "hbond")
cmd.set("dash_length", 0.20, "hbond")
cmd.set("dash_gap", 0.27, "hbond")

half_height = 9.15
rotation = camera(center, long_axis, vertical, half_height, 52, 23, roll=9, turn=-19)
width, height = 2350, 1200
ray("pocket.png", width, height)
selections = {
    "Trp286": "key and resi 286 and name CB",
    "Phe295": "key and resi 295 and name CZ",
    "Tyr337": "key and resi 337 and name CD2",
    "Trp86": "key and resi 86 and name CZ3",
    "donor": "ache and resi 295 and name N",
    "acceptor": "crystal and name O24",
}
anchors = {label: project(selection, center, rotation, half_height, width, height)
           for label, selection in selections.items()}
anchors.update(overview_ligand=overview_anchor, hbond_distance=hb_distance)
(OUT / "anchors.json").write_text(json.dumps(anchors, indent=2))
python end
