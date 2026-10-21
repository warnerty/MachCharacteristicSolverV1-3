import csv
from pathlib import Path
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from matplotlib.patches import Arc, Patch

# Developed by **Tyler Warner**

# read each ray and its flow state
path = Path(__file__).resolve().parent / "characteristics.csv"
lines = {}
with path.open(newline="") as file:
    for row in csv.DictReader(file):
        line = lines.setdefault(int(row["line_id"]), {
            "x": [], "y": [], "Mach": float(row["Mach"]),
            "theta": float(row["theta_deg"]), "alpha": float(row["alpha_deg"]),
            "pressure": float(row["pressure_Pa"]),
            "temperature": float(row["temperature_K"]),
        })
        line["x"].append(float(row["x"]))
        line["y"].append(float(row["y"]))

ordered = [lines[i] for i in sorted(lines)]
first, last = ordered[0], ordered[-1]
xc, yc = first["x"][0], first["y"][0]
theta = np.radians(last["theta"])
length = max(np.hypot(x - xc, y - yc)
             for line in ordered for x, y in zip(line["x"], line["y"]))
xmin, xmax = xc - 0.45 * length, xc + 1.08 * length
wall_end = yc + (xmax - xc) * np.tan(theta)
ymin = min(yc, wall_end) - 0.18 * length
ymax = max(max(line["y"]) for line in ordered) + 0.16 * length

# interpolate the pressure and temperature fields
X, Y = np.meshgrid(np.linspace(xmin, xmax, 650), np.linspace(ymin, ymax, 450))
angle = np.arctan2(Y - yc, X - xc)
order = np.argsort([line["alpha"] for line in ordered])
alpha = np.radians(np.array([line["alpha"] for line in ordered])[order])
pressure = np.array([line["pressure"] for line in ordered])[order] / 1000
temperature = np.array([line["temperature"] for line in ordered])[order]
P = np.interp(angle.ravel(), alpha, pressure).reshape(X.shape)
T = np.interp(angle.ravel(), alpha, temperature).reshape(X.shape)
wall_y = yc + np.maximum(X - xc, 0) * np.tan(theta)
mask = (Y < wall_y) | (np.hypot(X - xc, Y - yc) < 1e-10 * length)
P, T = np.ma.array(P, mask=mask), np.ma.array(T, mask=mask)

# draw the flow field and every characteristic
fig, ax = plt.subplots(figsize=(12, 7))
field = ax.pcolormesh(X, Y, P, cmap="viridis", shading="auto")
fig.colorbar(field, ax=ax, pad=0.03).set_label("Static pressure (kPa)")
for line in ordered:
    ax.plot(line["x"], line["y"], color="black", lw=0.6, alpha=0.35)

# label the temperature contours
if np.ptp(temperature) > 1e-8:
    levels = np.linspace(temperature.min(), temperature.max(), 9)[1:-1]
    contour = ax.contour(X, Y, T, levels=levels, colors="crimson",
                         linestyles="--", linewidths=1.3)
    labels = ax.clabel(contour, fmt=lambda value: f"{value:.1f} K", fontsize=8)
    for label in labels:
        label.set_bbox(dict(facecolor="white", edgecolor="none", alpha=0.8, pad=1))

# draw the wall and corner
ax.fill([xmin, xmin, xc, xmax, xmax], [ymin, yc, yc, wall_end, ymin],
        color="lightgray", zorder=4)
ax.plot([xmin, xc, xmax], [yc, yc, wall_end], color="black", lw=3, zorder=5)
ax.scatter(xc, yc, color="red", s=40, zorder=6)
ax.annotate("Expansion corner", xy=(xc, yc),
            xytext=(xc - 0.40 * length, yc - 0.12 * length),
            color="darkred", arrowprops=dict(arrowstyle="->", color="red"), zorder=7)
if abs(last["theta"]) > 0:
    ax.add_patch(Arc((xc, yc), 0.5 * length, 0.5 * length,
                     theta1=last["theta"], theta2=0, color="darkorange", lw=2, zorder=6))
    ax.text(xc + 0.35 * length, yc - 0.03 * length,
            f"Turn = {abs(last['theta']):g}°", color="darkorange", va="top", zorder=7)

# show inlet and outlet values inside the plot
box = dict(facecolor="white", edgecolor="none", alpha=0.85)
ax.text(xc - 0.35 * length, yc + 0.20 * length,
        f"Inlet: M = {first['Mach']:.3f}\n"
        f"{first['pressure'] / 1000:.1f} kPa, {first['temperature']:.1f} K",
        fontsize=9, bbox=box, zorder=7)
ax.text(xc + 0.55 * length, yc - 0.15 * length,
        f"Outlet: M = {last['Mach']:.3f}\n"
        f"{last['pressure'] / 1000:.1f} kPa, {last['temperature']:.1f} K",
        fontsize=9, bbox=box, zorder=7)
handles = [
    Line2D([], [], color="black", lw=3, label="Wall"),
    Line2D([], [], color="gray", lw=0.7, label="Mach characteristics"),
    Line2D([], [], color="crimson", ls="--", label="Static temperature contours (K)"),
    Patch(facecolor="lightgray", label="Solid region"),
]
ax.legend(handles=handles, loc="upper left", fontsize=8)
ax.set(xlabel="x (model units)", ylabel="y (model units)",
       title="Expansion Fan: Static Pressure and Temperature",
       xlim=(xmin, xmax), ylim=(ymin, ymax))
ax.set_aspect("equal", adjustable="box")
fig.tight_layout()
fig.savefig(path.parent / "pressure_temperature_field.png", dpi=180)
plt.show()
