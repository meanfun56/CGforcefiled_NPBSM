# Swarm-CG Optimization Data

This folder contains the input files, force-field parameters, and optimization outputs for the Swarm-CG parameterization of 24 natural product-derived small molecules within the Martini 3 framework.

## Directory structure
```text
Optimization_swarmCG/
├── README.md
├── BOR/
├── CAP/
├── CAT/
├── ...
└── (one folder per molecule, named by the molecule abbreviation)
```

Each molecule folder contains the files described below. For molecule-specific files, the prefix `XXX` denotes the molecule abbreviation (e.g., `BOR-md.tpr`, `CAP-md.tpr`).

## File description

| File | Description |
|------|-------------|
| `XXX-md.tpr` | GROMACS run input file for the all-atom molecular dynamics simulation. |
| `XXX.itp` | Molecule-specific topology/parameter file. |
| `aa_mapped_sasa.xvg` | Solvent-accessible surface area (SASA) data for the atomistic-mapped reference. |
| `cg_map.ndx` | Index file defining the atom-to-bead mapping. |
| `cg_model.itp` | Initial CG force-field parameters; also used as input for Swarm-CG optimization. |
| `cg_sasa.xvg` | SASA data for the CG model. |
| `equi.mdp` | GROMACS parameter file for equilibration. |
| `martini_v3.0.b.3.2.itp` | Martini 3.0 base force-field parameters. |
| `martini_v3.0_phospholipids.itp` | Martini 3.0 parameters for phospholipids (if applicable). |
| `martini_v3.0_solvents.itp` | Martini 3.0 solvent parameters. |
| `md.mdp` | GROMACS parameter file for production MD. |
| `mini.mdp` | GROMACS parameter file for energy minimization. |
| `opti_summary.png` | Summary of the EMD score evolution during iterative Swarm-CG optimization. |
| `optimized_CG_model_distributions.png` | Bonded distributions of the optimized CG model; the final EMD score is recorded here. |
| `reference_AA_distributions.png` | Reference bonded distributions from the atomistic simulations. |
| `start_conf.gro` | Initial CG system containing the CG molecule and water, used as input for Swarm-CG. |
| `system.top` | CG topology file required by Swarm-CG. |

## EMD scores

- The initial EMD score is recorded in the first distribution evaluation output (`distributions_eval_step_1.png` in the `all_evals_distributions` output).
- The final EMD score is recorded in `optimized_CG_model_distributions.png`.
- The full iterative optimization history is shown in `opti_summary.png`.

## Note

The all-atom trajectory file (`XXX-md-whole.xtc`) is not included because of its large file size (~800 MB per molecule). It is available from the corresponding author upon reasonable request.

## Related manuscript

These files support the parameterization and validation of the CG force field described in:

> Yan et al., *Coarse-Grained Force Field Parameters for Natural Product-Based Small Molecules: From Structural Mapping to Pharmaceutical Applications*.
