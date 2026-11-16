# Molecular Rotation: Laser-Induced Alignment of Linear Molecules

Python code that simulates the rotational dynamics of a linear molecule (for example CO₂) under a short, intense, non-resonant laser pulse. The molecule is treated quantum mechanically in the rigid-rotor basis, and the code outputs the time-dependent alignment parameter ⟨cos²θ⟩ before, during, and after the pulse.

## What it computes

- Thermal population of the rotational states *J*, weighted by the nuclear-spin statistics (`spin_even` / `spin_odd`).
- Time evolution of the rotational wavefunction under the interaction with a Gaussian laser pulse, using the polarizability anisotropy (`alpha_parallel`, `alpha_perpendicular`). The interaction term is `-¼ Δα E²(t)`, where `Δα = α∥ − α⊥`.
- Field-free rotational revival after the pulse, up to the chosen evaluation time.
- The thermal average of ⟨cos²θ⟩ as a function of time, written to a `.dat` file and plotted.

The time-dependent Schrödinger equation is integrated with `scipy.integrate.solve_ivp`. Angular matrix elements are computed with Wigner 3-j symbols from `sympy`.

## Repository layout

| Path | Contents |
|---|---|
| `Programs/` | Main working copy: the current solver, its input file, and example output. |
| `Programs/Rotation_linear.py` | **Current solver.** Reads `data.dat`, runs the simulation, writes the output file, and plots the result. |
| `Programs/data.dat` | Example input file for CO₂. |
| `Programs/Rotation_linear_old.py`, `..._v1.py`, `..._v2.py`, `v1.py`, `v2.py` | Earlier iterations of the solver. Kept for reference. |
| `Programs/solveivp.py`, `Programs/testpy.py`, `Programs/vren.py`, `Programs/test.m` | Small test and scratch scripts (ODE solver checks, a MATLAB `ode45` example). They are not part of the main calculation. |
| `Programs/*.dat`, `Programs/*.png` | Example outputs and plots from earlier runs. |
| `Linear-Molecular-Rotation/` | Separate git repository holding the same project with its own history. |

## Requirements

- Python 3.8 or newer
- [NumPy](https://numpy.org/)
- [SciPy](https://scipy.org/)
- [SymPy](https://www.sympy.org/)
- [Matplotlib](https://matplotlib.org/)

Install them with:

```bash
pip install numpy scipy sympy matplotlib
```

## Usage

1. Edit `data.dat` for your molecule and conditions (see the input format below).
2. Run the solver from the folder that contains `data.dat`. The script opens `data.dat` by relative path.

   ```bash
   cd Programs
   python Rotation_linear.py
   ```

3. The script prints progress and the convergence check on *J*, then:
   - writes `<Molecule>_at_<T>K_,<I0>_TWcm2.dat` in the working directory, and
   - opens a plot of ⟨cos²θ⟩ against time.

## Input file format (`data.dat`)

The solver reads one value per line. Only the first whitespace-, comma-, semicolon-, or hyphen-delimited token on each line is used, so anything after it is a free-form comment. Keep the order below.

| Line | Symbol in code | Meaning | Units |
|---|---|---|---|
| 1 | `Molecule_name` | Name used in output file names and the plot title | text |
| 2 | `B` | Rotational constant | cm⁻¹ |
| 3 | `alpha_par` | Polarizability parallel to the molecular axis | a.u. |
| 4 | `alpha_per` | Polarizability perpendicular to the molecular axis | a.u. |
| 5 | `spin_even` | Nuclear-spin degeneracy for even *J* | — |
| 6 | `spin_odd` | Nuclear-spin degeneracy for odd *J* (0 excludes odd *J*) | — |
| 7 | `temperature` | Rotational temperature | K |
| 8 | `J_MAX_user` | Maximum *J* included in the basis | — |
| 9 | `pulse_width` | Pulse duration (FWHM of intensity) | fs |
| 10 | `eval_time` | Total evaluation time after the pulse peak | fs |
| 11 | `wave_length` | Laser wavelength (read but not used in the current calculation) | nm |
| 12 | `I0` | Peak laser intensity | TW/cm² |

Example, the included CO₂ input:

```
CO2              = Molecule name
0.39021          =Rotational constant
27.25            =alpha_parallel
13.1             =alpha_perpendicular
1                =Spin_even
0                =Spin_odd
10               =Temperature
40               =J_max
10               =Pulse width in fs
300000            = Evaluation time in fs
800              =Wavelength
080              =Intensity in TW/cm^2
```

## Output

The `.dat` file has two whitespace-separated columns per row:

1. time (ps) on the solver's time grid
2. ⟨cos²θ⟩, thermally averaged and normalized

## Notes and limitations

- **Choose `J_MAX` large enough.** The script warns if population remains in the highest *J* state after the run. If it does, increase `J_MAX` in `data.dat`. It also reports the *J* value that captures 99.9 % of the thermal population.
- **Keep `spin_odd` consistent with the molecule.** For CO₂ (identical O nuclei, zero nuclear spin), odd *J* states are absent, so `spin_odd = 0`.
- **Intensity and pulse shape.** The field is Gaussian, with the pulse peak at `2 × pulse_width` in the internal time grid. The field is set to zero after `tend = 4 × pulse_width`.
- **Units.** Internal calculations use atomic units. Inputs are given in the units in the table above, and conversions are done in `standardize()`.
- Running the script opens a Matplotlib window. In a headless environment, comment out `molecule.plot()`.
- The `.dat` files in `Programs/` are outputs from earlier runs. Rerunning the solver writes new files with the same names and overwrites them.
