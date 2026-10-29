# Secondary-structure benchmark against DSSP

Each method is compared residue by residue with PyDSSP (3-state DSSP reimplementation, >97% agreement with the original DSSP according to its authors) on AlphaFold models of a stratified random sample of proteins from this app's dataset.

## app (DSSP port, proline-aware)

2026-10-24 · reference: PyDSSP 0.9.1 (numpy), 3-state DSSP reimplementation · 171 AlphaFold models (62,176 residues; 4 proteins without a model skipped)

| | Q3 agreement | helix F1 | strand F1 | loop F1 |
|---|---|---|---|---|
| Confident residues (pLDDT ≥ 70, n=48,650) | **99.6%** | 99.9% | 99.3% | 99.5% |
| All residues (n=62,176) | 99.7% | 99.9% | 99.2% | 99.7% |

Per-protein content correlation with DSSP: helix r = 1.0, strand r = 0.999

Confusion (confident residues; rows = DSSP, columns = app):

| DSSP \ app | H | E | C |
|---|---|---|---|
| H | 21,273 | 0 | 26 |
| E | 0 | 9,470 | 142 |
| C | 16 | 0 | 17,723 |

Q3 by protein class (confident residues): dna_binding 99.7%, ion_channel 99.9%, kinase 99.5%, other 99.6%, protease 99.7%, ribosomal_protein 99.3%, toxin 99.2%

## old C-alpha method (P-SEA-like)

2026-10-24 · reference: PyDSSP 0.9.1 (numpy), 3-state DSSP reimplementation · 171 AlphaFold models (62,176 residues; 4 proteins without a model skipped)

| | Q3 agreement | helix F1 | strand F1 | loop F1 |
|---|---|---|---|---|
| Confident residues (pLDDT ≥ 70, n=48,650) | **77.8%** | 86.3% | 74.9% | 66.9% |
| All residues (n=62,176) | 79.8% | 84.2% | 74.5% | 77.6% |

Per-protein content correlation with DSSP: helix r = 0.905, strand r = 0.96

Confusion (confident residues; rows = DSSP, columns = app):

| DSSP \ app | H | E | C |
|---|---|---|---|
| H | 19,884 | 49 | 1,366 |
| E | 602 | 7,655 | 1,355 |
| C | 4,320 | 3,131 | 10,288 |

Q3 by protein class (confident residues): dna_binding 82.2%, ion_channel 83.8%, kinase 75.7%, other 77.4%, protease 77.4%, ribosomal_protein 76.6%, toxin 81.3%
