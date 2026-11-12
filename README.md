# Vacuum-Propagation-for-STM-Simulation

`vpstm` provides two public workflows used in the manuscript:

1. **Q screening**: evaluates the near-source shape-and-vertical-derivative descriptor without doing a VP solve.
2. **VP**: propagates signed wavefunctions through a real-space total effective potential and returns wavefunction slices at requested tip heights.

The package deliberately uses portable NumPy arrays as the interface between CP2K and the analysis code. This keeps the solver independent of CP2K version-specific MOLog and Cube parsers.

## Install

```bash
python -m pip install -e .
```

## Q screening

Prepare `source_jet.npz` with arrays `F,Fz,Fzz,A,Az,Azz,weights,region,energies`. The first dimension of each field array is the sampled source-plane point and the second is the electronic state, i.e. shape `(n_point, n_state)`. `weights` and `energies` have shape `(n_state,)`; `region` is a nonnegative spatial weight normalized internally. `F*` are the DFT source-plane values and `A*` are the AVTA source-plane jet values evaluated on the same points. The `cp2k_outputs` block in the example records the CP2K provenance; it is not parsed implicitly.

```bash
vpstm examples/q_config.yaml
# equivalently: python -m vpstm examples/q_config.yaml
cat q_result.json
```

`Q` is a dimensionless near-source mismatch. It is a screening descriptor, not a VP image and not a replacement for a matched AVTA/VP validation on a new method or CP2K setup.

## VP propagation

Prepare (or export from the CP2K data-preparation step):

- `boundary.npy`: signed real or complex source wavefunction, shape `(ny,nx)`;
- `total_potential_eV.npy`: effective potential on the propagation grid, shape `(nz,ny,nx)`;
- `z_grid_A.npy`: the z coordinate of each potential slice.

```bash
vpstm examples/vp_config.yaml
```

The solver uses periodic lateral boundary conditions, a fixed source boundary, a zero Dirichlet upper boundary, and a second-order finite-difference discretization. It solves the stationary Schrödinger equation in the supplied potential and rejects an allowed-region interior unless `allow_allowed: true` is explicitly set. Complex source amplitudes are retained. The CLI writes both the complex propagated slices (`psi`) and their intensity maps (`intensity=|psi|^2`) to the output NPZ.

## Preparing CP2K data

The public interface expects a **signed wavefunction source boundary** and the **total local effective potential**. Standard CP2K Cube output may provide density, Hartree, XC, and local pseudopotential components separately, but it does not provide the total potential in the exact format required by this solver in every CP2K release.

This project used a small CP2K source modification to export the total potential during the SCF run. The relevant patches are in `cp2k_patch/`. In particular, `NATIVE_XC_CUBE.patch` adds the ordinary-GPW/PBE XC export used when constructing the total potential; its application and validation steps are summarized in `cp2k_patch/NATIVE_XC_CUBE.md`. Apply it only to the matching CP2K source revision, rebuild CP2K, and verify the output against the SCF energy and grid metadata before using it. Do not apply the patch blindly across CP2K versions: inspect the patch context and retain the patched executable hash in your provenance record.

A typical workflow is:

1. run a fixed-geometry CP2K SCF calculation;
2. export matched MO coefficients, the signed orbital/source data, density, and total potential using the patched executable;
3. interpolate the source boundary and total potential onto one common Cartesian grid;
4. run `vpstm` with `examples/vp_config.yaml`.

The repository does not claim that a single CP2K patch applies to every CP2K release. The patch is provided as a reproducibility aid for the revision used in this work. Because CP2K's MOLog and Cube formats vary across releases, the public solver intentionally separates CP2K export/interpolation from propagation; do not silently mix grids, spin channels, orbitals, or source planes.

## Data conventions

Energies are in eV, lengths in Å, and potentials in eV. The potential and boundary arrays must share the same lateral periodic cell and ordering. Complex wavefunctions are supported by the solver; Q source arrays may be complex. Keep the state window, spectral weights, source plane, lateral grid, and observation heights matched when comparing methods.

## Citation and license

This software accompanies the vacuum-propagation STM study. Add the manuscript citation after publication. See `LICENSE` for the code license.
