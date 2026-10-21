import csv
from pathlib import Path
import matplotlib.pyplot as plt

# Developed by **Tyler Warner**

# read coordinates beside this script
path = Path(__file__).resolve().parent / "characteristics.csv"
lines = {}
with path.open(newline="") as file:
    for row in csv.DictReader(file):
        line = lines.setdefault(int(row["line_id"]), {"x": [], "y": []})
        line["x"].append(float(row["x"]))
        line["y"].append(float(row["y"]))

# plot the characteristic rays
fig, ax = plt.subplots(figsize=(10, 6))
for line in lines.values():
    ax.plot(line["x"], line["y"], color="royalblue", linewidth=1)
ax.set(xlabel="x (model units)", ylabel="y (model units)",
       title="Prandtl–Meyer Expansion Fan")
ax.set_aspect("equal", adjustable="box")
ax.grid(True, alpha=0.3)
fig.tight_layout()
fig.savefig(path.parent / "expansion_fan.png", dpi=180)
plt.show()
