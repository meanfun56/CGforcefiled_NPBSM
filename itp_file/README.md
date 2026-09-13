# CG Force Field Parameters for 24 Natural Product-Derived Small Molecules

This folder contains the coarse-grained (CG) force-field parameter files (`.itp`) for the 24 natural product-derived small molecules parameterized in this work. These files are compatible with the Martini 3 force field and can be used directly in GROMACS simulations.

## File description

Each `.itp` file contains the bonded parameters (bond lengths, bond angles, dihedral angles, and force constants) for one molecule, as well as the bead type assignments. The non-bonded parameters for the bead types are taken from the standard Martini 3 force field.

The 24 molecules are:

| File | Molecule | File | Molecule |
|------|----------|------|----------|
| `FER.itp` | Ferulic acid | `BOR.itp` | Borneol |
| `EUG.itp` | Eugenol | `CNL.itp` | Eucalyptol |
| `CIN.itp` | Cinnamaldehyde | `GPD.itp` | Geniposide |
| `OST.itp` | Osthole | `CAP.itp` | Capsaicin |
| `IMP.itp` | Imperatorin | `RES.itp` | Resveratrol |
| `GEN.itp` | Genistein | `LIG.itp` | Ligustrazine |
| `PUE.itp` | Puerarin | `EPH.itp` | Ephedrine |
| `QUE.itp` | Quercetin | `MAG.itp` | Magnolol |
| `CAT.itp` | Catechin | `MET.itp` | Menthol |
| `EMO.itp` | Emodin | `PA.itp` | Perillaldehyde |
| `TPN.itp` | Terpinen-4-ol | `CNN.itp` | Limonene |
| `CIT.itp` | Citral | `LIN.itp` | Linalool |

## Usage

To use these parameters in a GROMACS simulation, include the corresponding `.itp` file in your system topology (`.top`) and ensure that the Martini 3 base force-field files (`martini_v3.0.b.3.2.itp`, `martini_v3.0_solvents.itp`, etc.) are also included. The bead types and non-bonded parameters follow the standard Martini 3 definitions.

Example:

```text
#include "martini_v3.0.b.3.2.itp"
#include "FER.itp"

[ molecules ]
FER 1
```
## Related manuscript
These parameter files support the findings in:

> Yan et al., Coarse-Grained Force Field Parameters for Natural Product-Based Small Molecules: From Structural Mapping to Pharmaceutical Applications.
