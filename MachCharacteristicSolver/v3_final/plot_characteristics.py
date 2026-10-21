import argparse
from pathlib import Path

import matplotlib.pyplot as plt
import matplotlib.patheffects as pe
import numpy as np
from matplotlib.collections import LineCollection
from matplotlib.colors import LinearSegmentedColormap
from matplotlib.lines import Line2D
from matplotlib.patches import Arc, Patch
from matplotlib.ticker import AutoMinorLocator, MaxNLocator

# Developed by **Tyler Warner**

def read_case(path):
    # read the current csv
    if not path.is_file():
        raise FileNotFoundError(f"Cannot find {path}. Put the CSV beside this script.")
    data = np.atleast_1d(np.genfromtxt(path, delimiter=",", names=True))
    required = {"line_id", "point_id", "x", "y", "Mach", "theta_deg",
                "alpha_deg", "pressure_Pa", "temperature_K"}
    missing = required - set(data.dtype.names or ())
    if missing:
        raise ValueError(f"Missing columns: {', '.join(sorted(missing))}. Regenerate the CSV.")
    origins = data[data["point_id"] == 0]
    if not len(origins) or np.max(data["point_id"]) <= 0:
        raise ValueError("The CSV needs corner points and ray endpoints.")
    corner = float(origins[0]["x"]), float(origins[0]["y"])
    rays = np.sort(data[data["point_id"] == np.max(data["point_id"])], order="line_id")
    if len(rays) < 2:
        raise ValueError("At least two rays are required.")
    for key in required:
        if not np.all(np.isfinite(rays[key])):
            raise ValueError(f"Nonfinite values in {key}.")
    if np.any(rays["pressure_Pa"] <= 0) or np.any(rays["temperature_K"] <= 0):
        raise ValueError("Pressure and temperature must be positive.")
    if not -90 < rays[-1]["theta_deg"] <= 0:
        raise ValueError("Use a downward lower-wall turn below 90 degrees.")
    return rays, corner


def endpoint(angle, bounds):
    # extend a straight ray to the frame
    xmin, xmax, ymin, ymax = bounds
    dx, dy = np.cos(angle), np.sin(angle)
    distances = []
    if dx > 1e-12:
        distances.append(xmax / dx)
    elif dx < -1e-12:
        distances.append(xmin / dx)
    if dy > 1e-12:
        distances.append(ymax / dy)
    elif dy < -1e-12:
        distances.append(ymin / dy)
    distance = min(s for s in distances if s > 0)
    return distance * dx, distance * dy


def draw_case(rays, corner, args):
    # colors and typography
    ink, muted, copper = "#17344d", "#667889", "#bd7339"
    solid_color = "#e9edf1"
    plt.rcParams.update({
        "font.family": "DejaVu Sans", "font.size": 9,
        "mathtext.fontset": "dejavusans", "axes.labelcolor": ink,
        "xtick.color": muted, "ytick.color": muted,
        "svg.fonttype": "none", "pdf.fonttype": 42, "savefig.facecolor": "white",
    })
    cmap = LinearSegmentedColormap.from_list(
        "static_pressure",
        ["#f5faf6", "#d3e9df", "#91cac2", "#4c9ca9", "#286987", "#123651"], N=512)

    # use coordinates relative to the corner
    rx, ry = rays["x"] - corner[0], rays["y"] - corner[1]
    length = float(np.max(np.hypot(rx, ry)))
    if not np.isfinite(length) or length <= 0:
        raise ValueError("The ray length must be positive and finite.")
    first, last = rays[0], rays[-1]
    theta = np.radians(float(last["theta_deg"]))
    turn = abs(float(last["theta_deg"]))
    xmin, xmax = -0.12 * length, 1.02 * length
    wall_end = xmax * np.tan(theta)
    ymin = min(0, wall_end) - 0.045 * length
    ymax = max(0.30 * length, float(np.max(ry))) + 0.045 * length
    bounds = xmin, xmax, ymin, ymax

    # interpolate pressure across the fan
    order = np.argsort(rays["alpha_deg"])
    alpha = np.radians(rays["alpha_deg"][order])
    pressure = rays["pressure_Pa"][order] / 1000
    temperature = rays["temperature_K"][order]
    X, Y = np.meshgrid(np.linspace(xmin, xmax, 1100), np.linspace(ymin, ymax, 850))
    angles = np.arctan2(Y, X)
    P = np.interp(angles.ravel(), alpha, pressure).reshape(X.shape)
    mask = (Y < np.maximum(X, 0) * np.tan(theta)) | (np.hypot(X, Y) < 1e-10 * length)
    P = np.ma.array(P, mask=mask)

    # size the page around the field
    aspect = (xmax - xmin) / (ymax - ymin)
    pw = min(7.4, 7.1 * aspect)
    ph = pw / aspect
    left, bottom = 0.72, 0.95
    fw, fh = left + pw + 3.0, bottom + ph + 1.65
    fig = plt.figure(figsize=(fw, fh), facecolor="white")

    def pos(x, y, w, h):
        return [x / fw, y / fh, w / fw, h / fh]

    ax = fig.add_axes(pos(left, bottom, pw, ph))
    cax = fig.add_axes(pos(left + pw + 0.20, bottom, 0.13, ph))
    side = fig.add_axes(pos(left + pw + 1.02, bottom, 1.85, ph))
    side.set_axis_off()
    side.set(xlim=(0, 1), ylim=(0, 1))
    ax.set_facecolor(solid_color)

    # title and subtitle
    fig.text(left / fw, (fh - 0.33) / fh, "GAS DYNAMICS  /  CENTERED EXPANSION FAN",
             color=muted, fontsize=8, weight="bold", va="top")
    fig.text(left / fw, (fh - 0.63) / fh, "Prandtl–Meyer expansion",
             color=ink, fontsize=23, weight="bold", va="top")
    fig.text(left / fw, (fh - 1.06) / fh,
             "Static pressure field with temperature isotherms", color=muted, fontsize=10, va="top")
    fig.text((fw - 0.20) / fw, (fh - 0.65) / fh, f"{turn:g}° TURN",
             color=copper, fontsize=15, weight="bold", ha="right", va="top")

    # pressure color scale
    field = ax.pcolormesh(X, Y, P, cmap=cmap, shading="nearest",
                         rasterized=True, vmin=0, vmax=float(np.max(pressure)), zorder=1)
    cbar = fig.colorbar(field, cax=cax)
    cbar.outline.set_visible(False)
    cbar.ax.tick_params(length=0, labelsize=8, pad=5)
    cbar.locator = MaxNLocator(nbins=6)
    cbar.update_ticks()
    cbar.ax.set_title("$p$\n(kPa)", color=ink, fontsize=9, pad=12)

    # show a limited selection of characteristics
    visible = np.unique(np.linspace(0, len(rays) - 1, min(args.rays, len(rays)), dtype=int))
    segments = [[(0, 0), endpoint(np.radians(rays[i]["alpha_deg"]), bounds)] for i in visible]
    ax.add_collection(LineCollection(segments, colors=ink, linewidths=0.45, alpha=0.19, zorder=2))

    # interpolate the local flow direction
    if args.streamlines:
        flow_angles = np.radians(rays["theta_deg"][order])
        theta_field = np.interp(angles.ravel(), alpha, flow_angles).reshape(X.shape)
        U = np.ma.array(np.cos(theta_field), mask=mask)
        V = np.ma.array(np.sin(theta_field), mask=mask)

        # start the streamlines above the upstream wall
        seed_y = np.linspace(0.04 * length, ymax - 0.04 * length, args.streamlines)
        seeds = np.column_stack((
            np.full(args.streamlines, xmin + 0.004 * length), seed_y
        ))
        first_arrow = len(ax.patches)
        streams = ax.streamplot(
            X[0, :], Y[:, 0], U, V,
            start_points=seeds,
            integration_direction="forward",
            broken_streamlines=False,
            density=2.0,
            color="white",
            linewidth=0.85,
            arrowsize=0.9,
            minlength=0.03,
            maxlength=10.0,
            zorder=2.5,
        )
        halo = [
            pe.Stroke(linewidth=1.7, foreground=ink, alpha=0.55),
            pe.Normal(),
        ]
        streams.lines.set_path_effects(halo)
        for arrow in ax.patches[first_arrow:]:
            arrow.set_path_effects(halo)

    # isotherms are straight rays in a centered fan
    if np.ptp(temperature) > 1e-8 and args.contours:
        levels = np.linspace(temperature.min(), temperature.max(), args.contours + 2)[1:-1]
        t_order = np.argsort(temperature)
        iso_angles = np.interp(levels, temperature[t_order], alpha[t_order])
        for i, (value, angle) in enumerate(zip(levels, iso_angles)):
            end = endpoint(angle, bounds)
            line, = ax.plot([0.018 * length * np.cos(angle), end[0]],
                            [0.018 * length * np.sin(angle), end[1]],
                            color=copper, lw=0.8, linestyle=(0, (5, 4)), zorder=3)
            line.set_path_effects([pe.Stroke(linewidth=1.65, foreground="white", alpha=0.55),
                                   pe.Normal()])
            radius = (0.73 if i % 2 == 0 else 0.91) * length
            ax.text(radius * np.cos(angle), radius * np.sin(angle), f"{value:.0f} K",
                    color="#874c25", fontsize=7.8, rotation=np.degrees(angle),
                    rotation_mode="anchor", ha="center", va="center", zorder=8,
                    bbox=dict(boxstyle="round,pad=0.18", facecolor="white", edgecolor="none", alpha=0.93))

    # wall and solid region
    ax.fill_between([xmin, 0, xmax], [ymin, ymin, ymin], [0, 0, wall_end],
                    color=solid_color, linewidth=0, zorder=4)
    wall_start = xmin + 0.012 * (xmax - xmin)
    ax.plot([wall_start, 0, xmax], [0, 0, wall_end],
            color=ink, lw=2, solid_capstyle="butt", zorder=5)
    ax.scatter(0, 0, s=28, facecolor=copper, edgecolor="white", linewidth=0.8, zorder=8)
    label_y = -(0.10 + 0.18 * abs(np.tan(theta))) * length

    # wall turning angle
    if turn:
        ax.plot([0, 0.23 * length], [0, 0], color="white", alpha=0.6,
                lw=0.7, linestyle=(0, (2, 3)), zorder=6)
        ax.add_patch(Arc((0, 0), 0.27 * length, 0.27 * length,
                         theta1=-turn, theta2=0, color=copper, lw=1.35, zorder=6))
        a = theta / 2
        ax.text(0.195 * length * np.cos(a), 0.195 * length * np.sin(a),
                rf"$\delta = {turn:g}^\circ$", fontsize=9, color="#874c25",
                va="center", ha="left", zorder=8,
                bbox=dict(boxstyle="round,pad=0.17", facecolor="white", edgecolor="none", alpha=0.9))

    # local flow directions
    ax.annotate("", xy=(0.04 * length, 0.15 * length),
                xytext=(xmin + 0.025 * length, 0.15 * length),
                arrowprops=dict(arrowstyle="-|>", color="white", lw=1.2, mutation_scale=11), zorder=7)
    ax.text(xmin + 0.025 * length, 0.19 * length, "UPSTREAM FLOW",
            fontsize=7.5, weight="bold", color="white", zorder=8)
    a = (theta + np.radians(float(last["alpha_deg"]))) / 2
    ox, oy = 0.66 * length * np.cos(a), 0.66 * length * np.sin(a)
    dx, dy = 0.095 * length * np.cos(theta), 0.095 * length * np.sin(theta)
    ax.annotate("", xy=(ox + dx, oy + dy), xytext=(ox - dx, oy - dy),
                arrowprops=dict(arrowstyle="-|>", color=ink, lw=1.25, mutation_scale=11), zorder=7)
    ax.annotate("EXPANDED FLOW", xy=(ox, oy), xytext=(0, 12), textcoords="offset points",
                color=ink, fontsize=7.5, ha="center", rotation=np.degrees(theta), zorder=8)

    # axes use the same physical scale
    unit = f" ({args.unit})" if args.unit else ""
    ax.set(xlim=(xmin, xmax), ylim=(ymin, ymax),
           xlabel="$x - x_c$" + unit, ylabel="$y - y_c$" + unit)
    ax.set_aspect("equal", adjustable="box")
    ax.xaxis.set_major_locator(MaxNLocator(nbins=6))
    ax.yaxis.set_major_locator(MaxNLocator(nbins=6))
    ax.xaxis.set_minor_locator(AutoMinorLocator(2))
    ax.yaxis.set_minor_locator(AutoMinorLocator(2))
    ax.tick_params(which="major", length=4, width=0.7, labelsize=8, pad=5)
    ax.tick_params(which="minor", length=2, width=0.5)
    ax.xaxis.labelpad, ax.yaxis.labelpad = 10, 10
    for name, spine in ax.spines.items():
        spine.set_color("#aab6c1" if name in ("top", "right") else "#738697")
        spine.set_linewidth(0.55)

    # keep numerical values outside the flow field
    def note(y, text, size=8, color=muted, weight="normal"):
        side.text(0, y, text, va="top", color=color, fontsize=size, weight=weight)

    note(0.995, "FLOW STATES", weight="bold")
    side.axhline(0.958, color="#dce3e9", lw=0.7)
    note(0.918, "UPSTREAM", weight="bold")
    note(0.866, f"$M_1$  {first['Mach']:.3f}", 18, ink, "bold")
    note(0.786, f"$p_1$   {first['pressure_Pa'] / 1000:,.1f} kPa", 10, ink)
    note(0.736, f"$T_1$   {first['temperature_K']:,.1f} K", 10, ink)
    side.axhline(0.674, color="#dce3e9", lw=0.7)
    note(0.634, "DOWNSTREAM", weight="bold")
    note(0.582, f"$M_2$  {last['Mach']:.3f}", 18, ink, "bold")
    note(0.502, f"$p_2$   {last['pressure_Pa'] / 1000:,.1f} kPa", 10, ink)
    note(0.452, f"$T_2$   {last['temperature_K']:,.1f} K", 10, ink)
    side.axhline(0.390, color="#dce3e9", lw=0.7)
    note(0.350, "CHANGE THROUGH FAN", weight="bold")
    changes = [("Mach", last["Mach"] / first["Mach"] - 1),
               ("Pressure", last["pressure_Pa"] / first["pressure_Pa"] - 1),
               ("Temperature", last["temperature_K"] / first["temperature_K"] - 1)]
    for y0, (name, value) in zip([0.294, 0.246, 0.198], changes):
        note(y0, name, 9)
        side.text(0.98, y0, f"{value * 100:+.1f}%", color=ink,
                  fontsize=10, weight="bold", ha="right", va="top")
    side.axhline(0.126, color="#dce3e9", lw=0.7)
    note(0.086, "PLANAR / ISENTROPIC", 7.5, weight="bold")
    note(0.044, r"Steady, inviscid, constant $\gamma$", 7.5)

    # compact legend above the field
    handles = [
        Line2D([], [], color=ink, lw=2, label="Wall"),
        Line2D([], [], color=muted, lw=0.7, label="Characteristics"),
        Line2D([], [], color=copper, lw=0.9, ls="--", label="Temperature (K)"),
        Patch(facecolor=solid_color, edgecolor="#c9d2dc", lw=0.5, label="Solid"),
    ]
    if args.streamlines:
        handles.insert(2, Line2D(
            [], [], color="white", lw=0.9, label="Streamlines",
            path_effects=[pe.Stroke(linewidth=1.7, foreground=ink), pe.Normal()],
        ))
    fig.legend(handles=handles, loc="lower left", frameon=False,
               bbox_to_anchor=(left / fw, (bottom + ph + 0.06) / fh),
               ncol=5, fontsize=8, handlelength=2.0, columnspacing=1.5, borderaxespad=0)

    # show the data count and units
    fig.text(left / fw, 0.36 / fh,
             f"{len(rays):,} computed rays  /  {len(visible)} displayed  "
             f"·  Ray length L = {length:g} {args.unit}".rstrip(), color=muted, fontsize=7.5)
    fig.text(left / fw, 0.17 / fh,
             "p: absolute static pressure   ·   T: static temperature   ·   "
             "Field interpolated between computed rays", color=muted, fontsize=7)
    return fig


def main():
    # edit the defaults or use command line options
    parser = argparse.ArgumentParser(description="Plot C++ expansion fan results.")
    parser.add_argument("csv", nargs="?",
                        default=str(Path(__file__).resolve().parent / "characteristics.csv"))
    parser.add_argument("--contours", type=int, default=10)
    parser.add_argument("--rays", type=int, default=40)
    parser.add_argument("--streamlines", type=int, default=7)
    parser.add_argument("--unit", default="model units")
    parser.add_argument("--no-show", action="store_true")
    parser.add_argument("--output", help="Output name without an extension.")
    args = parser.parse_args()
    if args.contours < 0 or args.rays < 2 or args.streamlines < 0:
        parser.error("Use at least two visible rays and nonnegative contour and streamline counts.")
    path = Path(args.csv).resolve()
    rays, corner = read_case(path)
    fig = draw_case(rays, corner, args)
    stem = Path(args.output).resolve() if args.output else path.parent / "pressure_temperature_field"
    stem.parent.mkdir(parents=True, exist_ok=True)

    # save a high resolution image and vector copies
    for extension in ("png", "pdf", "svg"):
        output = stem.with_suffix("." + extension)
        fig.savefig(output, dpi=450 if extension == "png" else 220)
        print(f"Saved: {output}")
    print(f"Reading: {path}")
    print(f"Turn: {abs(rays[-1]['theta_deg']):g} degrees | Inlet Mach: {rays[0]['Mach']:.3f}")
    if not args.no_show:
        plt.show()
    plt.close(fig)


if __name__ == "__main__":
    main()

