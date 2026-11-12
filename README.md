# v2qgt

[中文使用说明](README_Ch.md)

v2qgt calculates the quantum geometric tensor (QGT) from VASP wavefunctions. Based on [VASPBERRY](https://github.com/Infant83/VASPBERRY), it uses Kubo sums to calculate Berry curvature, the quantum metric, and a two-band QGT inequality.

The program reads the wavefunctions, energies, and k points from `WAVECAR`. The k-point count and order come from that file. If you change `KPOINTS`, generate a matching `WAVECAR` in VASP before running v2qgt.

## Build

The supplied Makefile uses Intel Fortran, Intel MPI, and MKL on Linux. Load the matching compiler and libraries, then run `make` in the source directory that contains the Makefile. Use `make clean && make` when changing compiler flags. Run the resulting `v2qgt.x` from the directory containing your `WAVECAR`, or pass the input path with `-f`.

## `-kubo` modes

| Mode | Quantity calculated |
| --- | --- |
| `0` (default) | Skips Kubo calculations. The legacy overlap-based Berry curvature calculation is disabled in the current main program. |
| `1` | Calculates band-resolved quantities in Cartesian `x,y`. Use `-qm 0` for Berry curvature or `-qm 1` for the quantum metric. A band range also produces a separate total. |
| `2`, `3`, `4` | Use the same calculation as mode `1`. |
| `5` | Evaluates the two-band QGT inequality using Cartesian `kx,ky` derivatives. Outputs the metric, traceless Berry curvature matrix, and both sides of the inequality at each k point. |
| `6` | Evaluates the same inequality using derivatives along two reciprocal lattice vectors selected by `QGT_KDIR1` and `QGT_KDIR2`. |
| `7` | Evaluates the same inequality using two arbitrary directions specified by `QGT_DVEC1` and `QGT_DVEC2`. |

Other positive integers currently fall through to the mode `1` calculation. Use `1` for that calculation.

## Command-line options

Pass options as flag–value pairs. Band numbers are one-based and follow the order in `WAVECAR`.

| Option | Modes 1–4 | Modes 5–7 |
| --- | --- | --- |
| `-nn N` | Band selection comes from `-is` or `-ii/-if`. | First band in the target pair. |
| `-mm M` | Second band for the `-abelian 2` calculation. | Second band in the target pair. Use a valid band number distinct from `N`. |
| `-qm 0/1` | `0` (default) outputs Berry curvature, with `xy` as the relevant component. `1` outputs `g_xx`, `g_xy`, and `g_yy`. | The routine computes both metric and curvature regardless of this flag. |
| `-abelian 1/2` | `1` (default) selects the single-band Abelian calculation. `2` selects a non-Abelian matrix element in a two-band subspace and requires `-mm`. | The routine always uses the two-band formulation. Specify `2` to acknowledge that choice and suppress the mode warning. |
| `-Qnm 0/1` | With `-abelian 2`, the intended selection is `0` for diagonal `Q_nn` or `1` for off-diagonal `Q_nm` in band space. See the source issue below. | The routine computes the full 2×2 band-space matrix internally. |
| `-is N` | Outputs only band `N`. | Select the pair with `-nn/-mm`. |
| `-ii N -if M` | Outputs each band from `N` through `M`, inclusive, and a separate total when the range contains multiple bands. | Select the pair with `-nn/-mm`. |

The modes 1–4 total uses a separate routine. It sums over the selected target bands with intermediate states above `-if`, so adding the individual band files can give a different result. The total routine is independent of the band-resolved `-abelian`, `-mm`, and `-Qnm` settings.

**Known source issue with `-Qnm`.** The main program declares `Q_nm` as an integer, while the parser and band-resolved non-Abelian routine implicitly declare their corresponding arguments as double precision. Make those types consistent and verify the result before relying on the modes 1–4 non-Abelian output.

| Option | Meaning |
| --- | --- |
| `-f FILE` | Input wavefunction file. Default is `WAVECAR`. |
| `-o TAG` | Output tag. For standard QGT runs, the file prefix is `BERRYCURV.TAG` when supplied and `BERRYCURV_KUBO` by default. |
| `-s 2` / `-s 1` | Read spinor (default, e.g. SOC) or non-spinor wavefunctions. Match the VASP calculation that produced the input. |
| `-kp N` | Controls Brillouin zone replication in modes 1–4. Default is `2`. Modes 5–7 write the original k points. |

## Directions and intermediate bands in modes 5–7

- In mode 6, set `QGT_KDIR1` and `QGT_KDIR2` to distinct integers from 1 through 3. They select from `b1`, `b2`, and `b3`. The defaults are 1 and 2.
- In mode 7, each direction is a quoted triple of coefficients in the reciprocal basis. For example, `QGT_DVEC1="1 1 0"` means `b1+b2`. The defaults are `"1 0 0"` and `"0 1 0"`. The two vectors must be nonzero and nonparallel. Their lengths affect the tensor components because the code uses the vectors as supplied. You can choose derivative directions independently of the sampling path in `WAVECAR`.
- `QGT_EXT_MODE=all` (default) includes every intermediate band outside the target pair. `near` includes only bands inside the inclusive window `QGT_EXT_BMIN` through `QGT_EXT_BMAX`. `remote` includes only bands outside that window. All three settings exclude the target bands.
- The Kubo denominator uses the mean energy of the two target bands as its reference. The program skips terms with a gap of `1e-6 eV` or less from that reference. Check `min_gap_eV` and `nskip_gap` for small denominators among the selected intermediate bands.

## Output files

Modes 1–4 write one file per selected band with `.EIG-N` in its name. The suffix is typically `.dat` for SOC or `.UP.dat` and `.DN.dat` for a collinear calculation with two spin channels. Each row contains three Cartesian k coordinates, the three result columns `xx,xy,yy`, and three fractional k coordinates. The header retains legacy Berry-curvature labels. Interpret the result columns according to `-qm`. The Chern number in this header has not been verified in this branch.

| Mode | Output suffix |
| --- | --- |
| `5` | `.qgtineq2d_kxky_MODE.dat` |
| `6` | `.qgtineq2d_k1k2_MODE.dat`, with the numbers replaced by the selected axes |
| `7` | `.qgtineq2d_dvec_MODE.dat` |

`MODE` is `all`, `near`, or `remote`. Use different `-o` tags for runs with different band pairs or mode 7 directions in the same directory.

Modes 5–7 write one row per input k point. The columns contain the k-point index and fractional coordinates, the half-trace metric components `g11,g22,g12`, `detg`, `lhs`, `omega_norm2`, `rhs`, `margin`, the real and imaginary parts of the four entries of the traceless Berry curvature matrix in order `11,12,21,22`, and `min_gap_eV,nskip_gap`.

The curvature convention is `Omega = i(Q12-Q21)`, with `Omega_SU2 = Omega - Tr(Omega) I/2`. The output uses `omega_norm2 = Tr(Omega_SU2²)/2`, `lhs = sqrt(max(detg,0))`, and `rhs = sqrt(max(omega_norm2,0))/2`. A nonnegative `margin = lhs - rhs` satisfies the inequality as evaluated by the code.

## Examples

Set `v2qgt` to the executable path and run from the directory with the input wavefunctions.

```bash
v2qgt=/path/to/v2qgt.x

# Band 10 Abelian Berry curvature and quantum metric
"$v2qgt" -kubo 1 -is 10 -qm 0 -o bc10
"$v2qgt" -kubo 1 -is 10 -qm 1 -o metric10

# Bands 10 and 11, with Cartesian kx,ky derivatives
"$v2qgt" -kubo 5 -nn 10 -mm 11 -abelian 2 -o pair10_11_xy

# Derivatives along b2,b3, with only bands 9 and 12 as intermediate states
QGT_KDIR1=2 QGT_KDIR2=3 \
QGT_EXT_MODE=near QGT_EXT_BMIN=9 QGT_EXT_BMAX=12 \
"$v2qgt" -kubo 6 -nn 10 -mm 11 -abelian 2 -o pair10_11_k2k3

# Derivatives along b1+b2 and b3
QGT_DVEC1="1 1 0" QGT_DVEC2="0 0 1" \
"$v2qgt" -kubo 7 -nn 10 -mm 11 -abelian 2 -o pair10_11_dvec
```

## MPI

Build from clean objects with `-DMPI_USE` in the Makefile's `F90` flags. Use the matching Intel MPI launcher. For example, four MPI ranks divide the k points among themselves and rank 0 writes the output file.

```bash
OMP_NUM_THREADS=1 MKL_NUM_THREADS=1 \
mpirun -np 4 "$v2qgt" -kubo 7 -nn 10 -mm 11 -abelian 2
```

## 📖 Citation
If you use **v2qgt** in your research, please cite our paper:

> **Sunje Kim**, **Yoonah Chung**, **Yuting Qian**, *et al.*,\
> **Direct measurement of the quantum metric tensor in solids**,\
> *Science*, **388**, 1050–1054 (2025).\
> [DOI: 10.1126/science.ado6049](https://www.science.org/doi/abs/10.1126/science.ado6049)

### BibTeX
```bibtex
@article{doi:10.1126/science.ado6049,
  author = {Sunje Kim and Yoonah Chung and Yuting Qian and Soobin Park and Chris Jozwiak and Eli Rotenberg and Aaron Bostwick and Keun Su Kim and Bohm-Jung Yang},
  title = {Direct measurement of the quantum metric tensor in solids},
  journal = {Science},
  volume = {388},
  number = {6751},
  pages = {1050-1054},
  year = {2025},
  doi = {10.1126/science.ado6049},
  url = {https://www.science.org/doi/abs/10.1126/science.ado6049},
  eprint = {https://www.science.org/doi/pdf/10.1126/science.ado6049}
}

@software{yutqian_2025_14927897,
  author       = {yutqian},
  title        = {yutqian/v2qgt: v2qgt v1.0.0-alpha.1},
  month        = feb,
  year         = 2025,
  publisher    = {Zenodo},
  version      = {v1.0.0-alpha.1},
  doi          = {10.5281/zenodo.14927897},
  url          = {https://doi.org/10.5281/zenodo.14927897}
}

```
