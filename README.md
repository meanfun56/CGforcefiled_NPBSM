# CG Force Field Parameters for Natural Product-Based Small Molecules

This repository contains the supporting data for the manuscript:

> **Coarse-Grained Force Field Parameters for Natural Product-Based Small Molecules: From Structural Mapping to Pharmaceutical Applications**  
> Yumei Yan, Xiaowen Wu, Minfang Feng, Xingxing Dai, Bingxuan Li, Xinhui Peng, Guoqing Li, Yueyang Xie, Zhixun Li, Rui Li, Xinyuan Shi*

## Repository structure

```text
CGforcefield_NPBSM/
├── README.md
├── AA2CG_mapping_rules/
│   ├── 24mol/
│   ├── figure1.png
│   ├── figure2.png
│   └── figure3.png
├── Optimization_swarmCG/
│   ├── BOR/
│   ├── CAP/
│   ├── CAT/
│   ├── ...
│   └── MAG/
├── PermeationMD/
│   ├── OST-CNN/
│   ├── OST-CNL/
│   ├── OST-pure/
│   ├── OST-TPN/
│   └── mdp/
├── b-CD_inclusionMD/
│   ├── FER-b-CD/
│   └── OST-b-CD/
├── script/
│   ├── gmx_interaction_rerun.py
│   └── gmx_membrane_analysis_v4.py
└── itp_file/
    ├── FER.itp
    ├── EUG.itp
    ├── ...
    └── MAG.itp
```

## Folder descriptions

### `AA2CG_mapping_rules/`
CG mapping schemes (bead mapping rules) for 24 natural product-derived small molecules.
- `24mol/`: individual mapping schemes for each molecule (e.g., `FER.png`, `EUG.png`).
- `figure1.png`–`figure3.png`: overview figures corresponding to Figures 1–3 in the main text.

### `Optimization_swarmCG/`
Swarm-CG optimization inputs and outputs for each molecule. One folder per molecule (e.g., `BOR/`, `CAP/`).
Each molecule folder contains:
- `cg_map.ndx`: atom-to-bead mapping index file.
- `cg_model.itp`: initial CG force-field parameters (Swarm-CG input).
- `XXX-md.tpr`: GROMACS run input file for the all-atom MD simulation.
- `XXX-md-whole.xtc`: corrected all-atom trajectory (not included due to large file size; available upon request).
- `mini.mdp`, `equi.mdp`, `md.mdp`: GROMACS parameter files for energy minimization, equilibration, and production MD.
- `start_conf.gro`: initial CG system containing the CG molecule and water (Swarm-CG input).
- `system.top`: CG topology file required by Swarm-CG.
- `martini_v3.0*.itp`: Martini 3.0 base force-field, phospholipid, and solvent parameters.
- `optimized_CG_model_distributions.png`: bonded distributions of the optimized CG model (final EMD score recorded here).
- `opti_summary.png`: EMD score evolution during iterative Swarm-CG optimization.
- `reference_AA_distributions.png`: reference bonded distributions from atomistic simulations.
- `MODEL_OPTI__STARTED_.../optimized_CG_model/`: optimized CG model output, containing `XXX.itp` (optimized CG parameters), `aa_mapped_sasa.xvg`, and `cg_sasa.xvg`.
- `model.../scg.log`: Rg and SASA values extracted from the optimization.

### `PermeationMD/`
Membrane permeation simulations. Each system (`OST-CNN/`, `OST-CNL/`, `OST-pure/`, `OST-TPN/`) contains three independent replicates (`1/`, `2/`, `3/`).
Each replicate contains:
- `step7_production.edr`, `step7_production.gro`, `step7_production.log`, `step7_production.tpr`.
The `mdp/` folder contains the GROMACS parameter files (`.mdp`) and associated output files for the sequential equilibration and production protocol.

### `b-CD_inclusionMD/`
β-cyclodextrin inclusion simulations. Each system (`FER-b-CD/`, `OST-b-CD/`) contains three independent replicates (`1/`, `2/`, `3/`).
Each replicate contains:
- `md.gro`, `md.log`, `md.tpr`, and `whole.xtc` (partial trajectory; complete trajectory available upon request).

### `script/`
Analysis scripts:
- `gmx_interaction_rerun.py`: re-runs GROMACS interaction energy calculations from existing trajectories.
- `gmx_membrane_analysis_v4.py`: analyzes membrane properties (APL, order parameters, diffusion coefficients, density distributions).

### `itp_file/`
Individual `.itp` files for the 24 molecules, copied from the optimized CG models. These files can be used directly in GROMACS simulations with the Martini 3 force field.

## Data availability

The datasets generated during and/or analysed during the current study are available in this GitHub repository:  
https://github.com/meanfun56/CGforcefield_NPBSM/tree/main

## Citation

If you use these data, please cite:

> Yan et al., *Coarse-Grained Force Field Parameters for Natural Product-Based Small Molecules: From Structural Mapping to Pharmaceutical Applications*.

## Contact

Xinyuan Shi  
School of Chinese Materia Medica, Beijing University of Chinese Medicine  
Email: xys_2019@126.com
```
