# β-CD Inclusion Molecular Dynamics Simulations

This folder contains the molecular dynamics (MD) simulation data for the inclusion of two natural product-derived drugs, ferulic acid (FER) and osthole (OST), into β-cyclodextrin (β-CD).

## Directory structure

```text
b-CD_inclusionMD/
├── README.md
├── FER-b-CD/
│   ├── 1/
│   │   ├── md.gro
│   │   ├── md.log
│   │   ├── md.tpr
│   │   └── whole.xtc
│   ├── 2/
│   │   └── ...
│   └── 3/
│       └── ...
└── OST-b-CD/
    ├── 1/
    │   └── ...
    ├── 2/
    │   └── ...
    └── 3/
        └── ...
```
Each drug–β-CD system (FER-b-CD and OST-b-CD) contains three independent replicate simulations, numbered 1, 2, and 3.

## File description
Each replicate folder contains the following GROMACS output files:

| File | Description |
|------|-------------|
| `md.gro` | Final configuration of the simulation. |
| `md.log` | GROMACS log file containing energy, temperature, pressure, and other runtime information. |
| `md.tpr` | GROMACS run input file (topology + parameters) used for the simulation. |
| `whole.xtc` | Partial trajectory used for analysis. |
## Note
md.gro — Final structure file in GROMACS .gro format, containing the coordinates of all atoms (or CG beads) at the end of the simulation, as well as the box dimensions. It can be used to visualize the final configuration.

md.log — GROMACS log file recording runtime information such as energies, temperature, pressure, density, and any warnings or errors. Useful for checking simulation stability and convergence.

md.tpr — Portable binary run input file for GROMACS. It contains the complete description of the simulation system (topology, force-field parameters, simulation parameters, and starting coordinates) and is required by gmx mdrun to run or extend the simulation.

whole.xtc — Partial trajectory in GROMACS compressed .xtc format, used for the analysis reported in the manuscript (e.g., calculating inclusion efficiency, radial distribution functions, and diffusion coefficients). The complete original trajectory is not included because of its large file size (~200 MB per replicate); it is available from the first author upon reasonable request.

These simulations were performed to investigate the inclusion behavior of FER and OST with β-CD and to calculate the inclusion efficiency. The results from the three replicates were averaged to obtain the mean ± standard deviation reported in Table 5 of the main manuscript.

## Related manuscript
These data support the findings in:

Yan et al., Coarse-Grained Force Field Parameters for Natural Product-Based Small Molecules: From Structural Mapping to Pharmaceutical Applications.
