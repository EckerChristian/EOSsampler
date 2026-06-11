# EOS/TOV Inference Code

This project samples dense-matter EOS tables, writes them to HDF5 files, and then solves the Tolman-Oppenheimer-Volkoff (TOV) equations for each accepted EOS. The code is organized as two MPI-aware executables:

- `EOS`: samples EOS tables, writes EOS-side parameters, and optionally writes an auxiliary pre-transition EOS extension.
- `TOV`: reads the EOS tables, solves stellar sequences, computes TOV-side likelihoods, and appends the results to the same HDF5 files.

Each MPI rank writes or processes its own HDF5 file:

```text
rank 0 -> 0.h5
rank 1 -> 1.h5
rank 2 -> 2.h5
...
```

Run the EOS and TOV stages with the same number of MPI ranks unless only a subset of rank files should be processed.

## Project layout

```text
EOScode/
  CMakeLists.txt
  README.md
  config/
    default.yaml
  include/
    constants.hpp
    config.hpp
    schema.hpp
    eos.hpp
    tov.hpp
    likelihoods.hpp
    hdf5_io.hpp
  src/
    main_eos.cpp
    main_tov.cpp
    config.cpp
    eos_sampler.cpp
    tov_solver.cpp
    likelihoods.cpp
    hdf5_io.cpp
  scripts/
    plot_results.py
```

External EOS input tables and likelihood files are not stored directly in the source tree. Their paths are configured in `config/default.yaml`.

## Dependencies

Required:

- C++17 compiler
- CMake 3.16 or newer
- MPI
- HDF5 with C API

Optional for plotting and inspection:

- Python 3
- `h5py`
- `numpy`
- `matplotlib`

On a cluster, load the corresponding modules before configuring the project, for example:

```bash
module load cmake
module load mpi
module load hdf5
```

## Configuration

The runtime configuration is stored in a flat YAML file. The default file is:

```text
config/default.yaml
```

The parser expects simple `key: value` entries. Do not use nested YAML blocks unless the parser is extended.

Typical development settings are:

```yaml
n_eos: 5
verbose: true

enable_mass_likelihood: true
enable_gw_likelihood: false
enable_xray_likelihoods: false
```

Typical production settings are:

```yaml
n_eos: 10000
verbose: false

enable_mass_likelihood: true
enable_gw_likelihood: true
enable_xray_likelihoods: true
```

The EOS stage uses the crust, CEFT, and pQCD input paths. The TOV likelihood stage additionally uses GW, NICER, and X-ray likelihood files. Check these paths before running:

```yaml
crust_eos_file: "../crust/BPS2.txt"
ceft_file: "../ceft/ceft2NLO_1.1ns.dat"
pqcd_file: "../pQCD/cs2_muB_mean_std_scale_ave.dat"
gw170817_file: "../constraints/LV/LV_prior.h5"
```

## Phase transitions and EOSext

If phase transitions are enabled, the EOS sampler stores the usual EOS table with the first-order transition included in the group:

```text
/<eos_index>/EOS
```

An optional auxiliary table can also be written to:

```text
/<eos_index>/EOSext
```

This table is controlled by:

```yaml
include_phase_transition: true
write_eos_extension: true
eos_extension_max_delta_mu_mev: 500.0
```

`EOSext` is constructed only by the `EOS` executable. The `TOV` executable reads the ordinary `/EOS` table and does not create or modify `/EOSext`.

The `EOSext` table contains the already constructed low-density EOS up to the left side of the phase transition. It then continues the last pre-transition sound-speed segment to larger chemical potential. The continuation stops at the first of the following conditions:

- the extrapolated sound speed reaches `cs2 = 1`,
- the extrapolated sound speed would reach `cs2 = 0`,
- the extension reaches `eos_extension_max_delta_mu_mev`.

The default maximum extension length is `500 MeV`. If the sound speed reaches `cs2 = 1`, the endpoint is finite and may be written. If the sound speed would reach `cs2 = 0`, the code stops before appending a singular point, because the density evolution contains the factor `1 / cs2`.

`EOSext` is intended as an auxiliary diagnostic/output table. It is not used by the TOV solver unless the TOV code is explicitly modified to read it.

## Building

From the project root:

```bash
rm -rf build
cmake -S . -B build -DCMAKE_BUILD_TYPE=Release
cmake --build build -j
```

This produces:

```text
build/EOS
build/TOV
```

For debugging, use:

```bash
rm -rf build
cmake -S . -B build -DCMAKE_BUILD_TYPE=Debug
cmake --build build -j
```

Release mode is recommended for production runs.

## Running

Run the EOS stage first:

```bash
cd build
mpirun -np 1 ./EOS ../config/default.yaml 12345
```

The final argument is the random seed. This creates `0.h5` for rank 0.

For multiple ranks:

```bash
mpirun -np 2 ./EOS ../config/default.yaml 12345
```

This creates:

```text
0.h5
1.h5
```

Then run the TOV stage with the same number of ranks:

```bash
mpirun -np 2 ./TOV ../config/default.yaml
```

The TOV executable opens the existing HDF5 files, reads each `/EOS` table, solves the TOV sequence, evaluates likelihoods, and writes the results back into the same files.

## HDF5 output schema

Each EOS sample is stored as a numbered group:

```text
/0
/1
/2
...
```

Each group contains:

```text
/<eos_index>/EOS
/<eos_index>/TOV
/<eos_index>/params
```

If `write_eos_extension: true` and a phase transition is present, the group also contains:

```text
/<eos_index>/EOSext
```

The ordinary EOS datasets are:

```text
/<eos_index>/EOS/e
/<eos_index>/EOS/p
/<eos_index>/EOS/n
/<eos_index>/EOS/mu
/<eos_index>/EOS/cs2
```

The optional EOS extension datasets are:

```text
/<eos_index>/EOSext/e
/<eos_index>/EOSext/p
/<eos_index>/EOSext/n
/<eos_index>/EOSext/mu
/<eos_index>/EOSext/cs2
```

The TOV datasets are:

```text
/<eos_index>/TOV/M
/<eos_index>/TOV/R
/<eos_index>/TOV/Lambda
/<eos_index>/TOV/stab
/<eos_index>/TOV/e_cent
/<eos_index>/TOV/mu_cent
/<eos_index>/TOV/n_cent
/<eos_index>/TOV/cs2_cent
/<eos_index>/TOV/Rqm
```

The parameter datasets include:

```text
/<eos_index>/params/dn
/<eos_index>/params/ePTl
/<eos_index>/params/pPT
/<eos_index>/params/nPTl
/<eos_index>/params/muPT
/<eos_index>/params/cs2PTr
/<eos_index>/params/cs2PTl
/<eos_index>/params/pCET
/<eos_index>/params/pPQCD
/<eos_index>/params/pM
/<eos_index>/params/pGW
/<eos_index>/params/ptot
/<eos_index>/params/nbranches
/<eos_index>/params/pXray
```

## Inspecting output

Use `h5ls`:

```bash
h5ls -r 0.h5 | head -100
```

Check whether `EOSext` was written:

```bash
h5ls -r 0.h5 | grep EOSext
```

Check whether TOV ran successfully:

```bash
h5ls -r 0.h5 | grep TOV
```

If `/0/TOV` exists but contains no datasets, the TOV stage likely failed before writing results. The EOS stage may create the empty `TOV` subgroup in advance.

Inspect the EOS extension numerically:

```bash
python3 - <<'PY'
import h5py

with h5py.File("0.h5", "r") as f:
    mu = f["0/EOSext/mu"][:]
    cs2 = f["0/EOSext/cs2"][:]
    print("EOSext length:", len(mu))
    print("mu range [MeV]:", mu[0], "->", mu[-1])
    print("cs2 range:", cs2.min(), "->", cs2.max())
    print("last points:")
    for m, c in zip(mu[-5:], cs2[-5:]):
        print(m, c)
PY
```

## Plotting results

If `scripts/plot_results.py` is present, run:

```bash
python ../scripts/plot_results.py 0.h5 --eos 0 --out plots
```

This creates plots such as:

```text
plots/0_eos0_p_of_e.png
plots/0_eos0_p_of_n.png
plots/0_eos0_cs2_of_mu.png
plots/0_eos0_mr.png
plots/0_eos0_lambda_of_m.png
```

Install Python plotting dependencies if needed:

```bash
pip install h5py numpy matplotlib
```

## Recommended validation workflow

Start with a small run:

```yaml
n_eos: 5
verbose: true
enable_mass_likelihood: true
enable_gw_likelihood: false
enable_xray_likelihoods: false
```

Then run:

```bash
cd build
mpirun -np 1 ./EOS ../config/default.yaml 12345
mpirun -np 1 ./TOV ../config/default.yaml
```

Inspect the output:

```bash
h5ls -r 0.h5 | head -100
h5ls -r 0.h5 | grep EOSext
h5ls -r 0.h5 | grep TOV
```

Then enable likelihoods in stages:

1. mass likelihood only
2. mass + GW likelihood
3. mass + GW + X-ray likelihoods

This makes it easier to identify file-path or grid-orientation issues.

## Common issues

### MPI abort on rank 1

If rank 1 fails immediately, check whether `1.h5` exists. The TOV executable expects one HDF5 file per MPI rank. If EOS was run with one rank, run TOV with one rank too:

```bash
mpirun -np 1 ./TOV ../config/default.yaml
```

If TOV should run with two ranks, first generate two EOS files:

```bash
mpirun -np 2 ./EOS ../config/default.yaml 12345
mpirun -np 2 ./TOV ../config/default.yaml
```

### Empty TOV group

An empty `TOV` subgroup usually means TOV failed before writing datasets. Rebuild in `Debug` mode or keep `verbose: true` and check the rank-specific error message.

### Missing EOSext group

If `/EOSext` is missing, check:

```yaml
include_phase_transition: true
write_eos_extension: true
eos_extension_max_delta_mu_mev: 500.0
```

Also make sure the `EOS` executable was rerun after changing the code or config. The `TOV` executable does not create `/EOSext`.

### Short EOSext continuation

A short continuation is not necessarily an error. The extension stops before the configured maximum if the extrapolated sound speed reaches `cs2 = 1` or would reach `cs2 = 0`.

### Duplicate pressure values

EOSs with first-order phase transitions may contain duplicate pressure values with different energy densities. This is expected. The TOV interpolation must allow non-decreasing pressure, not strictly increasing pressure.

### Non-increasing energy density

Energy density must be strictly increasing for the TOV-side interpolation. If TOV reports:

```text
EOS energy-density table must be strictly increasing for TOV interpolation
```

then the EOS sampler should reject that candidate before writing it. Regenerate the HDF5 files after fixing EOS validation.

### Likelihoods are exactly zero

If `pGW`, `pXray`, or `ptot` are exactly zero for all EOSs, common causes are:

- wrong likelihood file path
- wrong HDF5 dataset path
- radius and mass axes swapped
- lambda axes swapped
- TOV curve lies outside the likelihood grid

Test likelihood modules one at a time.

## Development notes

The intended dependency direction is:

```text
main_eos.cpp -> eos_sampler.cpp -> hdf5_io.cpp
main_tov.cpp -> tov_solver.cpp -> likelihoods.cpp -> hdf5_io.cpp
```

The EOS sampler should not depend on the TOV solver. The TOV solver should read only the HDF5 schema produced by the EOS sampler.

Physical constants belong in `include/constants.hpp`.

Runtime settings and file paths belong in `config/default.yaml` and `include/config.hpp`.

HDF5 group and dataset names belong in `include/schema.hpp`.

## License

Add a `LICENSE` file before distributing the code. For academic permissive reuse, BSD-3-Clause is a reasonable default, unless institutional or collaborator requirements specify another license.