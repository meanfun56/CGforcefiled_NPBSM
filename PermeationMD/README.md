# Permeation Molecular Dynamics Simulations

This folder contains the molecular dynamics (MD) simulation data for the transdermal permeation enhancement study. Four systems were simulated: pure membrane (OST-pure), and membranes containing three different permeation enhancers (OST-CNN, OST-CNL, OST-TPN). Each system was run in three independent replicates (1, 2, and 3).

## Directory structure

```text
PermeationMD/
├── README.md
├── OST-CNN/
│   ├── 1/
│   │   ├── step7_production.edr
│   │   ├── step7_production.gro
│   │   ├── step7_production.log
│   │   └── step7_production.tpr
│   ├── 2/
│   │   └── ...
│   └── 3/
│       └── ...
├── OST-CNL/
│   ├── 1/
│   │   └── ...
│   ├── 2/
│   │   └── ...
│   └── 3/
│       └── ...
├── OST-pure/
│   ├── 1/
│   │   └── ...
│   ├── 2/
│   │   └── ...
│   └── 3/
│       └── ...
├── OST-TPN/
│   ├── 1/
│   │   └── ...
│   ├── 2/
│   │   └── ...
│   └── 3/
│       └── ...
└── mdp/
    ├── step6.0_minimization.edr
    ├── step6.0_minimization.gro
    ├── step6.0_minimization.log
    ├── step6.0_minimization.mdp
    ├── step6.0_minimization.tpr
    ├── step6.0_minimization.trr
    ├── step6.1_minimization.edr
    ├── step6.1_minimization.gro
    ├── step6.1_minimization.log
    ├── step6.1_minimization.mdp
    ├── step6.1_minimization.tpr
    ├── step6.1_minimization.trr
    ├── step6.2_equilibration.cpt
    ├── step6.2_equilibration.edr
    ├── step6.2_equilibration.gro
    ├── step6.2_equilibration.log
    ├── step6.2_equilibration.mdp
    ├── step6.2_equilibration.tpr
    ├── step6.2_equilibration.xtc
    ├── step6.4_equilibration.xtc
    ├── step6.5_equilibration.cpt
    ├── step6.5_equilibration.edr
    ├── step6.5_equilibration.gro
    ├── step6.5_equilibration.log
    ├── step6.5_equilibration.mdp
    ├── step6.5_equilibration.tpr
    ├── step6.5_equilibration.xtc
    ├── step6.6_equilibration.cpt
    ├── step6.6_equilibration.edr
    ├── step6.6_equilibration.gro
    ├── step6.6_equilibration.log
    ├── step6.6_equilibration.mdp
    ├── step6.6_equilibration.tpr
    ├── step6.6_equilibration.xtc
    ├── step7_production.edr
    ├── step7_production.gro
    ├── step7_production.log
    └── step7_production.tpr
```

## File description

Each `OST-XXX` system folder contains three replicate runs (`1`, `2`, `3`). Each replicate contains the production simulation output files:

| File | Description |
|------|-------------|
| `step7_production.edr` | Energy file from the production run. |
| `step7_production.gro` | Final configuration from the production run. |
| `step7_production.log` | GROMACS log file from the production run. |
| `step7_production.tpr` | Run input file for the production run. |

The `mdp/` folder contains the GROMACS parameter files (`.mdp`) and associated output files for the sequential equilibration and production protocol:

| File | Description |
|------|-------------|
| `step6.0_minimization.*` | First energy minimization step. |
| `step6.1_minimization.*` | Second energy minimization step. |
| `step6.2_equilibration.*` | First equilibration step. |
| `step6.4_equilibration.xtc` | Trajectory from an intermediate equilibration step. |
| `step6.5_equilibration.*` | Equilibration step. |
| `step6.6_equilibration.*` | Final equilibration step. |
| `step7_production.*` | Production run (reference files). |

File extensions follow standard GROMACS conventions: `.mdp` (parameter file), `.tpr` (run input), `.gro` (coordinates), `.edr` (energies), `.log` (log), `.cpt` (checkpoint), `.xtc` (compressed trajectory), `.trr` (full-precision trajectory).

## Note

The complete production trajectories (`step7_production.xtc`) are not included due to their large file size. They are available from the first author upon reasonable request.

## Related manuscript

These data support the findings in:

> Yan et al., *Coarse-Grained Force Field Parameters for Natural Product-Based Small Molecules: From Structural Mapping to Pharmaceutical Applications*.


