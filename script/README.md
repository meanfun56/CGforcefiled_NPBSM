# Analysis Scripts for CG-MD Simulations

This folder contains the Python scripts used for analyzing the coarse-grained molecular dynamics (CG-MD) simulation trajectories in this work.

## File description

| File | Description |
|------|-------------|
| `gmx_interaction_rerun.py` | Script for re-running GROMACS interaction energy calculations (e.g., van der Waals interaction energies between drug–drug, drug–β-CD, drug–membrane, and drug–solvent pairs) from existing trajectories. |
| `gmx_membrane_analysis_v4.py` | Script for analyzing membrane properties, including area per lipid (APL), lipid order parameters (S), diffusion coefficients (D), and density distributions along the membrane normal (z-axis). |

## Usage

Both scripts are written in Python and require a working GROMACS installation (version 2020.6 or later) and the `gmx` command-line tools. They also depend on standard Python libraries such as `numpy`, `matplotlib`, and `MDAnalysis` (or `gmx` output parsing utilities).

Typical usage examples:

```bash
# Calculate interaction energies
python gmx_interaction_rerun.py -s system.tpr -f trajectory.xtc -n index.ndx -o interaction_energy.xvg

# Analyze membrane properties
python gmx_membrane_analysis_v4.py -s system.tpr -f trajectory.xtc -n index.ndx -o membrane_analysis/
```
## Related manuscript
These data support the findings in:

> Yan et al., *Coarse-Grained Force Field Parameters for Natural Product-Based Small Molecules: From Structural Mapping to Pharmaceutical Applications*.
