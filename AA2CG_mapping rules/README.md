# CG Mapping Schemes (Bead Mapping Rules)

This folder contains the coarse-grained (CG) mapping schemes for the 24 natural product-derived small molecules parameterized in this work.

## Directory structure
```text
AA2CG_mapping_rules/
├── README.md
├── 24mol/
│   ├── FER.png
│   ├── EUG.png
│   ├── ...
│   └── MAG.png
├── figure1.png
├── figure2.png
└── figure3.png
```

The `24mol/` folder contains individual mapping scheme images for each of the 24 molecules, named by the molecule abbreviation. The three `figure*.png` files provide overviews of the mapping rules at different levels.

## File description

| File | Description |
|------|-------------|
| `24mol/` | Individual CG mapping schemes for each of the 24 molecules. Each image shows the all-atom structure and the corresponding Martini 3 CG beads. |
| `figure1.png` | Overview of the bead mapping rules for different functional groups (corresponds to Figure 1 in the main text). |
| `figure2.png` | Overview of the CG mapping schemes for different small-molecule core scaffolds (corresponds to Figure 2 in the main text). |
| `figure3.png` | Overview of the CG mapping schemes for all 24 natural product-derived small molecules (corresponds to Figure 3 in the main text). |

## Note

- The 24 molecules are: FER, EUG, CIN, OST, IMP, GEN, PUE, QUE, CAT, EMO, MET, PA, TPN, CNN, CIT, LIN, BOR, CNL, GPD, CAP, RES, LIG, EPH, MAG.
- The mapping scheme for MAG is presented in Figure 4 of the main text.
- All mapping schemes follow the hierarchical construction rules described in the main text.

## Related manuscript

These data support the findings in:

> Yan et al., *Coarse-Grained Force Field Parameters for Natural Product-Based Small Molecules: From Structural Mapping to Pharmaceutical Applications*.
```
