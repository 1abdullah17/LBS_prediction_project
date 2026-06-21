from pymol import cmd
from PIL import Image, ImageDraw
import os

# =====================================================
# FILES
# =====================================================

CRYSTAL_FILE = "data/pdb_clean/5P9J_clean.pdb"
AF_FILE = "data/alphafold_clean/AF-Q06187-F1_clean.pdb"

OUTPUT_DIR = "figures"

os.makedirs(
    OUTPUT_DIR,
    exist_ok=True
)

# =====================================================
# IMAGE SETTINGS
# =====================================================

WIDTH = 2000
HEIGHT = 1500

# =====================================================
# COMMON VIEW
# =====================================================

cmd.reinitialize()

cmd.load(CRYSTAL_FILE, "crystal")
cmd.load(AF_FILE, "af")

cmd.align(
    "af and name CA",
    "crystal and name CA"
)

cmd.orient()

VIEW = cmd.get_view()

cmd.delete("all")

# =====================================================
# PANEL A
# Crystal + ligand
# =====================================================

cmd.load(CRYSTAL_FILE, "crystal")

cmd.hide("everything")

cmd.show("cartoon", "polymer")

cmd.color("grey80", "polymer")

cmd.select(
    "ligand",
    "resn 8E8"
)

cmd.show(
    "sticks",
    "ligand"
)

cmd.color(
    "red",
    "ligand"
)

cmd.set_view(VIEW)

cmd.bg_color("white")

cmd.set(
    "ray_opaque_background",
    0
)

cmd.ray(
    WIDTH,
    HEIGHT
)

panelA = f"{OUTPUT_DIR}/panelA.png"

cmd.png(
    panelA,
    dpi=300
)

cmd.delete("all")

# =====================================================
# PANEL B
# Crystal vs AlphaFold
# =====================================================

cmd.load(CRYSTAL_FILE, "crystal")
cmd.load(AF_FILE, "af")

rmsd = cmd.align(
    "af and name CA",
    "crystal and name CA"
)[0]

cmd.hide("everything")

cmd.show(
    "cartoon",
    "crystal"
)

cmd.show(
    "cartoon",
    "af"
)

cmd.color(
    "grey70",
    "crystal"
)

cmd.color(
    "cyan",
    "af"
)

cmd.set_view(VIEW)

cmd.bg_color("white")

cmd.set(
    "ray_opaque_background",
    0
)

cmd.ray(
    WIDTH,
    HEIGHT
)

panelB = f"{OUTPUT_DIR}/panelB.png"

cmd.png(
    panelB,
    dpi=300
)

cmd.delete("all")

# =====================================================
# PANEL C
# AlphaFold + transferred ligand
# =====================================================

cmd.load(AF_FILE, "af")
cmd.load(CRYSTAL_FILE, "crystal")

cmd.align(
    "crystal and name CA",
    "af and name CA"
)

cmd.create(
    "transferred_ligand",
    "crystal and resn 8E8"
)

cmd.hide("everything")

cmd.show(
    "cartoon",
    "af"
)

cmd.color(
    "cyan",
    "af"
)

cmd.show(
    "sticks",
    "transferred_ligand"
)

cmd.color(
    "red",
    "transferred_ligand"
)

cmd.set_view(VIEW)

cmd.bg_color("white")

cmd.set(
    "ray_opaque_background",
    0
)

cmd.ray(
    WIDTH,
    HEIGHT
)

panelC = f"{OUTPUT_DIR}/panelC.png"

cmd.png(
    panelC,
    dpi=300
)

# =====================================================
# COMBINE PANELS
# =====================================================

imgA = Image.open(panelA)
imgB = Image.open(panelB)
imgC = Image.open(panelC)

panel_w = 1400
panel_h = 1000

imgA = imgA.resize(
    (panel_w, panel_h)
)

imgB = imgB.resize(
    (panel_w, panel_h)
)

imgC = imgC.resize(
    (panel_w * 2, panel_h)
)

canvas = Image.new(
    "RGB",
    (panel_w * 2, panel_h * 2),
    "white"
)

canvas.paste(
    imgA,
    (0, 0)
)

canvas.paste(
    imgB,
    (panel_w, 0)
)

canvas.paste(
    imgC,
    (0, panel_h)
)

draw = ImageDraw.Draw(canvas)

draw.text(
    (30, 30),
    "A",
    fill="black"
)

draw.text(
    (panel_w + 30, 30),
    "B",
    fill="black"
)

draw.text(
    (30, panel_h + 30),
    "C",
    fill="black"
)

final_file = (
    f"{OUTPUT_DIR}/BTK_three_panel_figure.png"
)

canvas.save(
    final_file,
    dpi=(600, 600)
)

print("\n================================")
print("RMSD =", round(rmsd, 3), "Å")
print("Saved:", final_file)
print("================================")

cmd.quit()