# Molecular Dynamics Simulation of Multi-Particle Electrostatic Cold Spray Deposition of Aluminum Using EAM Potential

<p align="center">
  <img src="https://img.shields.io/badge/LAMMPS-MD%20Simulation-blue?style=for-the-badge&logo=gnu&logoColor=white"/>
  <img src="https://img.shields.io/badge/Aluminum-EAM%20Potential-orange?style=for-the-badge"/>
  <img src="https://img.shields.io/badge/Cold%20Spray-Electrostatic-red?style=for-the-badge"/>
  <img src="https://img.shields.io/badge/Multi--Particle-Deposition-purple?style=for-the-badge"/>
  <img src="https://img.shields.io/badge/OVITO-VTK%20Export-green?style=for-the-badge"/>
  <img src="https://img.shields.io/badge/License-CC%20BY--NC%204.0-lightgrey?style=for-the-badge"/>
</p>

<p align="center">
  A fully atomistic <b>molecular dynamics simulation of multi-particle electrostatic cold spray deposition</b>
  of aluminum using LAMMPS. Models two spherical FCC Al nanoparticles accelerated by a uniform
  electrostatic field toward an Al substrate, capturing impact dynamics, stress wave propagation,
  and plastic deformation via the <i>EAM/alloy potential (Al99.eam.alloy)</i>. Von Mises stress
  and hydrostatic pressure are computed per-atom across approach, impact, and post-impact relaxation phases.
</p>
<img width="1600" height="1200" alt="multi ecs" src="https://github.com/user-attachments/assets/0959af6b-04b9-4209-ac75-a819efcabf1f" />

---

## Physics

The simulation captures the full sequence of electrostatic cold spray deposition at the atomistic scale:

- FCC aluminum substrate and spherical particles modeled with the EAM/alloy interatomic potential
- Electrostatic field acceleration of charged nanoparticles toward the substrate
- Per-atom charge assignment and field-force coupling via `fix efield`
- NVE integration for particles and NVT (Nose-Hoover) thermostat for the substrate
- Per-atom stress tensor computation via `compute stress/atom`
- Von Mises stress and hydrostatic pressure fields extracted across all phases
- Multi-particle simultaneous release to model sequential impact and particle-particle interaction

---

## Geometry

```
  Simulation box (units: Angstrom):

  z = 110  +------------------------------------------+  <- box top
           |                                          |
  z =  73  |        [ Particle 2 ]  center (0,0,63)  |  <- radius 10 Ang
           |                                          |
  z =  50  |        [ Particle 1 ]  center (0,0,40)  |  <- radius 10 Ang, gap ~3 Ang
           |                                          |
  z =   5  +==========================================+  <- substrate top
           |            SUBSTRATE                    |
  z = -15  +------------------------------------------+  <- frozen bottom layer
  z = -20  +------------------------------------------+  <- box bottom

  Box XY: -50 to +50 Angstrom (periodic in x and y)
  Box Z:    -20 to +110 Angstrom (shrink-wrapped in z)
```

- **Substrate**: FCC Al slab, z = -20 to +5 Ang; bottom 5 Ang frozen (zero force)
- **Particle 1**: Sphere, center (0, 0, 40), radius = 10 Ang
- **Particle 2**: Sphere, center (0, 0, 63), radius = 10 Ang (3 Ang gap above Particle 1)
- **Boundary**: Periodic in x, y; shrink-wrapped (s) in z

---

## Simulation Phases

| Phase | Run Steps | Simulated Time | Description |
|-------|-----------|----------------|-------------|
| Equilibration | 1,000 | 1 ps | System thermalization before release |
| Approach | 2,000 | 2 ps | Both particles travel toward substrate |
| Impact | 4,000 | 4 ps | High-resolution impact capture (dump every 50 steps) |
| Post-impact relaxation | 4,000 | 4 ps | E-field removed; system relaxes |
| Final relaxation | 3,000 | 3 ps | Further equilibration after deposition |
| **Total** | **14,000** | **14 ps** | — |

---

## Material Parameters

### Aluminum Properties

| Parameter | Symbol | Value | Unit |
|-----------|--------|-------|------|
| Crystal structure | — | FCC | — |
| Lattice constant | a | 4.05 | Angstrom |
| Atomic mass | m | 26.9815 | g/mol |
| Interatomic potential | — | EAM/alloy | Al99.eam.alloy |

### Electrostatic Parameters

| Parameter | Symbol | Value | Unit |
|-----------|--------|-------|------|
| Applied electric field | Ez | -5.0e-3 | V/Angstrom |
| Particle charge per atom | q | 0.006 | e (elementary charge) |
| Substrate charge | — | 0.0 | e |
| Initial particle velocity | vz | -3.0 | Angstrom/ps |

### Thermostat Settings

| Group | Fix | Temperature | Damping |
|-------|-----|-------------|---------|
| Particles (both) | NVE | — | — |
| Substrate | NVT (Nose-Hoover) | 300 K | 0.1 ps |
| Bottom layer | Frozen | — | setforce 0 |

---

## Governing Equations

### Electrostatic Force on Charged Atoms

```
F_z = q * Ez

where:
  q   = charge per atom (e)
  Ez  = applied electric field component (V/Ang)
  F_z = force per atom in z-direction (eV/Ang)
```

### EAM/Alloy Potential Energy

```
E_i = F_alpha(sum_{j != i} rho_alpha_beta(r_ij)) + 0.5 * sum_{j != i} phi_alpha_beta(r_ij)

where:
  F_alpha    = embedding energy function
  rho        = electron density contribution
  phi        = pair interaction potential
```

### Per-Atom Stress Tensor (LAMMPS units: bar * Ang^3)

```
sigma_ab = (1/V_i) * [ -m_i * v_ia * v_ib + 0.5 * sum_j (r_iab * f_ijb) ]
```

### Von Mises Stress

```
sigma_VM = sqrt( 0.5 * [ (sxx-syy)^2 + (syy-szz)^2 + (szz-sxx)^2
                        + 6*(sxy^2 + sxz^2 + syz^2) ] )
```

### Hydrostatic Pressure

```
P = -(sxx + syy + szz) / 3
```

---

## Group Definitions

```
substrate   <- FCC Al slab atoms (z = -20 to +5 Ang)
particle    <- Particle 1 atoms  (sphere center (0,0,40), r=10)
particle2   <- Particle 2 atoms  (sphere center (0,0,63), r=10)
particles   <- union of particle + particle2  (used for efield and NVE)
bottom      <- frozen bottom layer (z = -20 to -15 Ang)
mobile      <- all atoms EXCEPT bottom (substrate mobile + both particles)
```

---

## Output Files

| File | Dump Frequency | Content |
|------|---------------|---------|
| `deposition_approach.lammpstrj` | every 500 steps | id, type, x, y, z, full stress tensor, von Mises, pressure |
| `vonmises.lammpstrj` | every 500 steps | id, type, x, y, z, von Mises only |
| `impact_detail.lammpstrj` | every 50 steps | id, type, x, y, z, von Mises, pressure (high-res impact) |
| `deposition_postimpact.lammpstrj` | every 500 steps | id, type, x, y, z, full stress tensor, von Mises, pressure |
| `final_deposited_ecs_2p.data` | — | Final atomic positions and velocities (LAMMPS data format) |
| `deposition_ecs_2p.restart` | — | Binary restart file for continuation runs |
| `frequencies.txt` | — | Final printed stress statistics |

### Console Output (Final Analysis Block)

```
================================================
=== ELECTROSTATIC COLD SPRAY RESULTS        ===
=== TWO-PARTICLE DEPOSITION                 ===
================================================
E-field applied (V/Ang):         -0.005
Particle charge per atom (e):    0.006
Maximum von Mises stress: <value> bars
Average von Mises stress: <value> bars
Maximum pressure (compression):  <value> bars
Minimum pressure (tension):      <value> bars
================================================
```

---

## Repository Structure

```
cold_spray_md/
|
|-- electrostatic_cold_spray_2particles.lammps   # Main LAMMPS simulation script
|-- Al99.eam.alloy                               # EAM potential file (required)
|-- README.md                                    # This file
|
|-- outputs/                                     # Generated on run
    |-- deposition_approach.lammpstrj            # Approach phase trajectory
    |-- vonmises.lammpstrj                       # Von Mises stress trajectory
    |-- impact_detail.lammpstrj                  # High-resolution impact frames
    |-- deposition_postimpact.lammpstrj          # Post-impact relaxation
    |-- final_deposited_ecs_2p.data              # Final atomic configuration
    |-- deposition_ecs_2p.restart                # Restart file
```

---

## How to Run

### Requirements

- LAMMPS (any recent version with EAM and KSPACE packages): https://www.lammps.org
- EAM potential file `Al99.eam.alloy` (available from NIST IPRP or LAMMPS potentials directory)
- OVITO or VMD for trajectory visualization: https://www.ovito.org

### Step 1 — Obtain the EAM potential file

```bash
# Download Al99.eam.alloy from the LAMMPS potentials directory or NIST IPRP
# Place it in the same directory as the script
ls Al99.eam.alloy
```

### Step 2 — Run the simulation

```bash
lmp -in electrostatic_cold_spray_2particles.lammps
```

Or in parallel with MPI:

```bash
mpirun -np 4 lmp -in electrostatic_cold_spray_2particles.lammps
```

The script will:
1. Build the FCC substrate and two spherical nanoparticles
2. Minimize energy to remove initial overlaps
3. Assign charges and initial velocities to both particles simultaneously
4. Run equilibration, approach, impact, and relaxation phases
5. Compute and print von Mises stress and pressure statistics
6. Write final data and restart files

### Step 3 — Visualize in OVITO

1. Open OVITO > `File > Load File` > select `deposition_approach.lammpstrj`
2. Add modifier: **Color Coding** → set property to `v_von_mises`
3. Add modifier: **Slice** to view cross-section of impact zone
4. Step through frames to observe particle approach and impact
5. Open `impact_detail.lammpstrj` for higher time-resolution playback of the impact event

---

## What to Look for in Results

### Impact Dynamics

Particle 1 impacts the substrate first. Particle 2 follows approximately 23 Ang behind. This models the effect of a leading particle pre-deforming the substrate before the trailing particle arrives — a key mechanism in sequential cold spray deposition.

### Von Mises Stress Distribution

Peak von Mises stress concentrates at the particle-substrate interface during impact. High-stress zones propagate as elastic waves through the substrate. Post-impact stress relaxes as kinetic energy dissipates into the thermostatted substrate.

### Pressure Field

Compressive pressure (positive values) dominates beneath the impact zone. Tensile regions (negative pressure) appear at the particle periphery and substrate surface, indicating potential void nucleation sites in real spray conditions.

### Particle Bonding

Successful deposition is indicated by persistent atom-atom contacts across the original particle-substrate interface in the final configuration. The `write_data` output can be re-loaded to measure adhesion geometry.

---

## Common Errors and Fixes

| Error | Cause | Fix |
|-------|-------|-----|
| `Cannot open potential file Al99.eam.alloy` | Potential file not found | Place `Al99.eam.alloy` in the run directory |
| `Atoms set too far from proc subdomain` | Box too small for particle positions | Extend z-upper limit in box region |
| `Lost atoms` during impact | Atoms ejected beyond simulation box | Add `thermo_modify lost warn` or extend z boundary |
| Double integration warning | NVE and NVT applied to same atoms | Ensure `fix 1` uses `particles` group and `fix 2` uses `substrate` group only |
| Dump file already exists error | Old dump files from previous run present | Delete or rename old `.lammpstrj` files before re-running |
| High initial temperature spike | Particle-substrate overlap at t=0 | Run `minimize` before dynamics (already included in script) |
| Particles not moving | E-field direction incorrect | Verify `Ez_field` is negative for downward acceleration |

---

## Extending the Model

| Extension | What to Change |
|-----------|----------------|
| Three or more particles | Add `particle3`, `particle4` regions; extend box z; union all into `particles` group |
| Different impact velocity | Change `velocity particle set 0.0 0.0 -X.X` |
| Offset particle positions | Change sphere center x, y coordinates (e.g., `sphere 10 0 40 10`) |
| Different material (Cu, Ni) | Change `mass`, `pair_coeff` potential file, and lattice constant |
| Higher temperature substrate | Change NVT target in `fix 2 substrate nvt temp T T 0.1` |
| Larger particles | Increase sphere radius; extend box z accordingly |
| Random particle placement | Use `create_atoms random` with a loop over particle regions |
| Radial distribution function | Add `compute rdf all rdf 100` and `fix ave/time` after deposition |

---

## Citation

If you use this code in your research, please cite:

```bibtex
@software{mishra_2026_coldspray,
  author    = {Mishra, A.},
  title     = {Molecular Dynamics Simulation of Multi-Particle Electrostatic Cold Spray Deposition of Aluminum Using EAM Potential},
  year      = {2026},
  publisher = {Zenodo},
  doi       = {10.5281/zenodo.20061969},
  url       = {https://doi.org/10.5281/zenodo.20061969}
}
```

Plain text citation:

> Mishra, A. (2026). *Molecular Dynamics Simulation of Multi-Particle Electrostatic Cold Spray Deposition of Aluminum Using EAM Potential*. Zenodo. https://doi.org/10.5281/zenodo.20061969

---

## Author

**akshansh11**
GitHub: https://github.com/akshansh11

---

## License

<p>
<a rel="license" href="http://creativecommons.org/licenses/by-nc/4.0/">
<img alt="Creative Commons Licence" style="border-width:0" src="https://i.creativecommons.org/l/by-nc/4.0/88x31.png"/>
</a>
<br/>
This work is licensed under a
<a rel="license" href="http://creativecommons.org/licenses/by-nc/4.0/">Creative Commons Attribution-NonCommercial 4.0 International License</a>.
</p>

You are free to:

- **Share** — copy and redistribute the material in any medium or format
- **Adapt** — remix, transform, and build upon the material

Under the following terms:

- **Attribution** — You must give appropriate credit to akshansh11 and provide a link to this repository
- **NonCommercial** — You may not use the material for commercial purposes

Copyright 2026 akshansh11. All rights reserved for commercial use.
