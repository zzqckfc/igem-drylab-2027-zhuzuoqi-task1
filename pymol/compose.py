"""Lay out the ray-traced PyMOL images, legends, and residue labels."""

import json
import os
from pathlib import Path

from PIL import Image, ImageDraw, ImageFont


ROOT = Path(__file__).resolve().parents[1]
RENDER = ROOT / "pymol" / "_render"
FIGURES = ROOT / "figures"
INK = "#172b39"
MUTED = "#435967"
RULE = "#d9e3e8"
TEAL = "#007a75"
ORANGE = "#e84f17"
PURPLE = "#8a468e"
WIDTH = 2760
MARGIN = 106
PANEL_X = (120, 1430)


def font(size, bold=False):
    windows = Path(os.environ.get("WINDIR", "C:/Windows")) / "Fonts"
    candidates = (
        windows / ("segoeuib.ttf" if bold else "segoeui.ttf"),
        "DejaVuSans-Bold.ttf" if bold else "DejaVuSans.ttf",
    )
    for candidate in candidates:
        try:
            return ImageFont.truetype(str(candidate), size)
        except OSError:
            continue
    return ImageFont.load_default(size)


def text(draw, position, value, size, color=INK, bold=False, anchor="la"):
    draw.text(position, value, font=font(size, bold), fill=color, anchor=anchor)


def plate(height, title, subtitle):
    image = Image.new("RGB", (WIDTH, height), "white")
    draw = ImageDraw.Draw(image)
    text(draw, (MARGIN, 57), title, 82, bold=True)
    text(draw, (MARGIN, 165), subtitle, 42, MUTED)
    rule(draw, 249)
    return image, draw


def rule(draw, y):
    draw.line((MARGIN, y, WIDTH - MARGIN, y), fill=RULE, width=3)


def render(image, name, xy):
    with Image.open(RENDER / f"{name}.png") as layer:
        molecule = layer.convert("RGBA")
    image.paste(molecule, xy, molecule)
    return molecule.getchannel("A")


def paired_views(image, draw, views, top, height):
    """Keep the frame dimensions, baselines and offsets equal across both views."""
    frames = []
    for x, tag, (caption, name) in zip(PANEL_X, ("a", "b"), views):
        text(draw, (x, 278), tag, 49, bold=True)
        text(draw, (x + 64, 284), caption, 43, bold=True)
        alpha = render(image, name, (x, top))
        if alpha.size != (1210, height):
            raise ValueError(f"Inconsistent panel frame for {name}")
        left, upper, right, lower = alpha.getbbox() or (0, 0, 0, 0)
        if min(left, upper, 1210-right, height-lower) <= 0:
            raise ValueError(f"Empty or cropped molecular render for {name}")
        frames.append(alpha)
    return frames


def swatch(draw, x, y, label, color):
    draw.rounded_rectangle((x, y + 17, x + 49, y + 31), radius=7, fill=color)
    text(draw, (x + 66, y), label, 40)


def save(image, filename):
    image.save(FIGURES / filename, dpi=(300, 300), optimize=True)


def shifted(point, offset):
    return tuple(round(value + delta) for value, delta in zip(point, offset))


def label_bounds(draw, name, xy, size, side, alpha, offset):
    bounds = draw.textbbox(xy, name, font=font(size, True), anchor=side)
    left, top, right, bottom = bounds
    dx, dy = offset
    overlap = alpha.crop((left-dx-10, top-dy-10, right-dx+10, bottom-dy+10))
    if overlap.getextrema()[1] > 20:
        raise ValueError(f"{name} label overlaps the molecular render")
    return bounds


def leader(draw, name, label_xy, atom_xy, *, alpha, offset, side="la", color=INK):
    bounds = label_bounds(draw, name, label_xy, 49, side, alpha, offset)
    x = bounds[2] + 17 if atom_xy[0] > bounds[2] else bounds[0] - 17
    y = (bounds[1] + bounds[3]) // 2
    line_color = MUTED if color == INK else color
    draw.line((x, y, *atom_xy), fill=line_color, width=4)
    tip_x, tip_y = atom_xy
    draw.ellipse((tip_x-5, tip_y-5, tip_x+5, tip_y+5), fill=line_color)
    text(draw, label_xy, name, 49, color, bold=True, anchor=side)


def overview():
    image, draw = plate(1685, "Donepezil in human acetylcholinesterase",
                        "4EY7 · Chain A · X-ray structure at 2.35 Å")
    top = 355
    alpha, _ = paired_views(image, draw, [("Full protein fold", "overview"),
                                         ("Local binding-gorge surface", "gorge")], top, 1130)
    anchors = json.loads((RENDER / "anchors.json").read_text())
    offset = (PANEL_X[0], top)
    tip = shifted(anchors["overview_ligand"], offset)
    leader(draw, "E20", (120, 1040), tip, color=TEAL, alpha=alpha, offset=offset)
    rule(draw, 1520)
    swatch(draw, 106, 1564, "Crystal donepezil (E20)", TEAL)
    swatch(draw, 810, 1564, "Residues within 6 Å of E20", "#547e9b")
    swatch(draw, 1774, 1564, "AChE molecular surface", "#a4bdcb")
    save(image, "figure1_overview.png")


def overlay():
    image, draw = plate(1320, "Crystal and re-docked donepezil",
                        "4EY7 · Both poses shown in the crystal receptor coordinate frame")
    paired_views(image, draw, [("Side view", "overlay_a"), ("90° about the ligand axis",
                                "overlay_b")], 355, 720)
    rule(draw, 1095)
    swatch(draw, 106, 1168, "Crystal pose", TEAL)
    swatch(draw, 586, 1168, "Vina rank 1", ORANGE)
    result = json.loads((ROOT / "docking" / "result.json").read_text())
    rmsd = result["crystal_heavy_atom_rmsd_angstrom"]
    right = WIDTH - MARGIN
    text(draw, (right, 1136), f"Heavy-atom RMSD  {rmsd:.3f} Å", 53, bold=True, anchor="ra")
    text(draw, (right, 1210), f"{result['heavy_atoms']} matched heavy atoms",
         38, MUTED, anchor="ra")
    save(image, "figure2_superposition.png")


def pocket():
    image, draw = plate(1660, "Contacts in the AChE binding gorge",
                        "4EY7 · Selected residues around crystallographic donepezil")
    offset = (205, 277)
    alpha = render(image, "pocket", offset)
    anchors = json.loads((RENDER / "anchors.json").read_text())

    # Perimeter labels stay clear of the molecular models.
    labels = {
        "Trp286": (550, 1200, "ra"),
        "Phe295": (1810, 345, "la"),
        "Tyr337": (2110, 560, "la"),
        "Trp86": (2190, 1155, "la"),
    }
    for name, (x, y, side) in labels.items():
        leader(draw, name, (x, y), shifted(anchors[name], offset),
               side=side, alpha=alpha, offset=offset)

    middle = [(start + end) / 2 for start, end in zip(anchors["donor"], anchors["acceptor"])]
    distance_xy = shifted(middle, (offset[0] - 53, offset[1] - 31))
    distance = f"{anchors['hbond_distance']:.2f} Å"
    label_bounds(draw, distance, distance_xy, 44, "ra", alpha, offset)
    text(draw, distance_xy, distance, 44, PURPLE, bold=True, anchor="ra")

    rule(draw, 1490)
    swatch(draw, 106, 1540, "Crystal donepezil (E20)", TEAL)
    for x in range(1390, 1450, 23):
        draw.rounded_rectangle((x, 1564, x + 14, 1570), radius=3, fill=PURPLE)
    text(draw, (1468, 1537), "Hydrogen bond · Phe295 N ··· E20 O24", 44)
    save(image, "figure3_binding_pocket.png")


if __name__ == "__main__":
    FIGURES.mkdir(parents=True, exist_ok=True)
    overview()
    overlay()
    pocket()
