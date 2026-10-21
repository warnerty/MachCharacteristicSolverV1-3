# Mach Characteristic Solver
Developed by **Tyler Warner**

A C++ and Python tool for calculating and visualizing a two-dimensional, centered Prandtl–Meyer expansion fan.

The solver calculates the flow state on each characteristic, exports coordinates and thermodynamic properties to CSV, and produces a figure with a static-pressure field, temperature isotherms, curved streamlines, and an inlet/downstream summary. The current C++ implementation uses OpenMP to calculate independent characteristic states.

![Current expansion fan visualization](docs/images/current_streamlines.png)

## What the project demonstrates

- Compressible-flow calculations using the Prandtl–Meyer relation.
- Numerical inversion with a bracketed bisection solver.
- Characteristic directions and coordinate generation.
- Isentropic static pressure and temperature.
- Parallel state calculation with OpenMP.
- Separation of numerical calculations in C++ from visualization in Python.
- Scientific figure design and PNG, PDF, and SVG export.

This implementation handles a single centered expansion fan. It is a foundation for further method-of-characteristics work; interacting waves, characteristic intersections, and reflections are future extensions.

## Physical model

The model assumes steady, planar, inviscid, isentropic flow of a calorically perfect gas with constant specific-heat ratio `gamma`.

Flow initially travels in the positive x direction. A lower wall turns clockwise through a positive expansion magnitude `delta`, and the fluid occupies the region above the wall.

For each turning increment:

$$
\theta_i=-\delta_i,\qquad
\nu_i=\nu(M_1)+\delta_i
$$

The Prandtl–Meyer function is

$$
\nu(M)=
\sqrt{\frac{\gamma+1}{\gamma-1}}
\tan^{-1}\sqrt{\frac{\gamma-1}{\gamma+1}(M^2-1)}
-\tan^{-1}\sqrt{M^2-1}.
$$

After solving for Mach number, the Mach angle and fan-ray direction are

$$
\mu_i=\sin^{-1}(1/M_i),\qquad
\alpha_i=\theta_i+\mu_i.
$$

Coordinates along a ray follow

$$
x=x_c+s\cos\alpha_i,\qquad
y=y_c+s\sin\alpha_i.
$$

Pressure and temperature are calculated relative to the inlet **static** conditions:

$$
\frac{T_i}{T_1}=
\frac{1+\frac{\gamma-1}{2}M_1^2}
{1+\frac{\gamma-1}{2}M_i^2},
\qquad
\frac{p_i}{p_1}=
\left(\frac{T_i}{T_1}\right)^{\gamma/(\gamma-1)}.
$$

Angles use radians internally and degrees in the CSV. Pressure inputs are absolute, in pascals; temperature inputs are in kelvin.

## Build and run

Requirements: a C++17 compiler with OpenMP support, Python, NumPy, and Matplotlib. GCC/MSYS2 is the intended Windows workflow.

Open a terminal **in the repository root** and install the plotting dependencies into your selected Python environment:

~~~powershell
python -m pip install -r requirements.txt
~~~

Compile:

~~~powershell
g++ -O3 -std=c++17 -fopenmp characteristicSolver.cpp -o characteristicSolver.exe
~~~

After compilation succeeds, run:

~~~powershell
.\characteristicSolver.exe
~~~

Generate the figure:

~~~powershell
python plot_characteristics.py
~~~

On Linux, compile to `characteristicSolver` and run `./characteristicSolver`.

C++ prints the absolute path of the CSV it writes. Its output directory is the terminal's current working directory. Python reads `characteristics.csv` beside the script by default. Run both programs from the repository root to keep these paths aligned.

## Configure a case

Edit the input constants in `characteristicSolver.cpp`, save, rebuild, and rerun C++ before plotting.

| Input | Packaged default | Meaning |
|---|---:|---|
| `gamma` | 1.4 | Constant ratio of specific heats |
| `inletMach` | 2.0 | Inlet Mach number |
| `turnDegrees` | 25.0 | Positive magnitude of the downward wall turn |
| `inletPressure` | 200000.0 | Absolute static pressure, Pa |
| `inletTemperature` | 300.0 | Static temperature, K |
| `cornerX`, `cornerY` | 0.0, 0.0 | Expansion-corner coordinates |
| `rayLength` | 10.0 | Length used to generate ray endpoints |
| `intervals` | 100000 | Turning intervals; creates 100001 rays |
| `segmentsPerRay` | 1 | One segment; two endpoints per ray |

The input Mach must be at least one, `gamma` must exceed one, and the requested Prandtl–Meyer angle must remain within its finite-Mach range. The plotting geometry supports a downward lower-wall turn below 90°.

Length units are consistent with the user's model. If the C++ coordinates are in meters, label them accordingly:

~~~powershell
python plot_characteristics.py --unit m
~~~

For a useful parameter comparison, keep the inlet conditions fixed and try turn angles of 5°, 15°, and 25°.

### Reproducing the packaged example

Running the packaged default case produced these downstream values:

| Quantity | Value |
|---|---:|
| Mach number | 3.0855 |
| Static pressure | 37.4944 kPa |
| Static temperature | 185.9469 K |

The bundled example images were generated from the packaged source. The basic and thermodynamic snapshots use the same inlet Mach and wall turn, so their flow states can be compared directly.

## Plot controls

~~~powershell
python plot_characteristics.py --streamlines 10 --contours 12 --rays 25
~~~

| Option | Default | Effect |
|---|---:|---|
| `--streamlines` | 7 | Number of seeded streamline paths; 0 hides them |
| `--contours` | 10 | Number of temperature isotherms; 0 hides them |
| `--rays` | 40 | Maximum number of visible characteristic rays |
| `--unit` | model units | Axis-unit label |
| `--no-show` | off | Export without opening the figure window |
| `--output` | pressure_temperature_field | Output filename stem |
| Positional CSV path | beside script | Read a different CSV file |

For example, a separate input and output location can be supplied:

~~~powershell
python plot_characteristics.py examples/run.csv --output figures/run --no-show
~~~

All computed states contribute to interpolation even when fewer rays are displayed. Pressure is represented by color, temperature by dashed isotherms, and streamlines by solid curves with arrows.

The plot coordinates are relative to the expansion corner: `x - x_c` and `y - y_c`. Straight characteristic rays are extended to the plotting frame. The color field and streamlines are reconstructed from the exported angular states.

## Streamlines

Python interpolates the local flow angle `theta_deg` over the plotting grid and integrates the direction field

$$u_{*}=\cos\theta,\qquad v_{*}=\sin\theta.$$

These unit direction components are sufficient to trace streamline geometry. They do not provide dimensional velocity magnitudes or particle travel times.

In a centered fan, Mach, pressure, and temperature are constant along a characteristic ray and vary across the fan. Streamlines follow the local velocity direction and turn as they cross the characteristic rays.

## Outputs

| Output | Contents |
|---|---|
| `characteristics.csv` | Ordered characteristic samples and flow states |
| `pressure_temperature_field.png` | 450 dpi figure |
| `pressure_temperature_field.pdf` | Figure with vector annotations and a rasterized pressure field |
| `pressure_temperature_field.svg` | Figure with editable vector elements and an embedded pressure field |

The current CSV columns are:

~~~text
line_id,point_id,x,y,Mach,theta_deg,nu_deg,mu_deg,alpha_deg,pressure_Pa,temperature_K
~~~

With the default inputs, the CSV contains 200002 data rows. Two endpoints fully describe each straight ray; exporting more intermediate points does not change its geometry.

## Parallelization and resolution

Only the independent characteristic-state calculation is parallelized. Each thread writes to a separate element of a preallocated array. CSV export remains sequential to preserve ordering and avoid concurrent writes to one file.

Increasing `intervals` refines the angular sampling. Increasing `rayLength` changes the spatial extent without increasing the number of states or samples. CSV row count is

$$
(N+1)(S+1),
$$

where `N` is the number of turning intervals and `S` is the number of segments per ray.

The Python field grid has its own resolution. More characteristic states alone do not increase that grid's dimensions. Runtime also depends on CSV writing, CSV loading, and figure rendering; parallel speedup is not guaranteed for the complete workflow.

GCC must receive `-fopenmp` for both the directives and OpenMP runtime functions. The program reports calculation and export times separately.

## Development progression

Two self-contained earlier stages are included:

| Stage | Folder | Main capability |
|---|---|---|
| v1: basic fan | [history/v1_basic_fan](history/v1_basic_fan) | Sequential characteristic generation and a simple ray plot |
| v2: thermodynamic fields | [history/v2_thermodynamics](history/v2_thermodynamics) | Static pressure/temperature, wall geometry, contours, and a legend |
| Current | Repository root | OpenMP, reduced CSV redundancy, refined figures, and streamlines |

**Snapshot provenance:** the earlier folders are reconstructed functional stages based on code shared during development. Their inputs were standardized to one comparison case, and minor packaging/input-check changes keep them runnable. They are not recovered Git commits or byte-exact copies of unshared local files.

[View the illustrated iteration history](docs/ITERATION_HISTORY.md).

### Earlier visualizations

| Basic characteristic fan | Thermodynamic fields |
|---|---|
| ![Basic fan](docs/images/v1_basic_fan.png) | ![Thermodynamic fields](docs/images/v2_thermodynamics.png) |

## Files

| Path | Purpose |
|---|---|
| `characteristicSolver.cpp` | Current parallel C++ solver |
| `plot_characteristics.py` | Current plotting and streamline script |
| `requirements.txt` | Python dependencies |
| `REPOSITORY_DESCRIPTION.txt` | Short description for the repository About field |
| `history/` | Runnable reconstructed earlier stages |
| `docs/ITERATION_HISTORY.md` | Feature comparison and progression figures |
| `docs/images/` | Reproducible example images used in this README |

## Scope and future work

The model treats one ideal expansion fan. It does not solve a general two-dimensional CFD mesh, shock waves, boundary layers, real-gas effects, or interactions between wave families.

Useful next steps include characteristic-intersection calculations, wall-reflection boundary conditions, independent reference-case comparisons, and measured scaling studies.

## Author

**Tyler Warner**  
Mechanical Engineering Student — University of Connecticut

Developed as a personal engineering project exploring Prandtl–Meyer expansion fans, numerical methods, OpenMP parallelization, and scientific visualization.

[GitHub](https://github.com/warnerty)

## References

- [NASA Glenn: Centered Expansion Fan](https://www.grc.nasa.gov/WWW/BGH/expans.html)
- [NASA Glenn: Isentropic Flow Equations](https://www.grc.nasa.gov/WWW/BGH/isentrop.html)
- [GCC: Enabling OpenMP](https://gcc.gnu.org/onlinedocs/libgomp/Enabling-OpenMP.html)
- [Matplotlib: streamplot](https://matplotlib.org/stable/api/_as_gen/matplotlib.pyplot.streamplot.html)


