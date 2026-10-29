# Protein function classifier - evaluation report

Generated 2026-10-24T22:55:19. All data downloaded from the **live UniProtKB REST API**.

## Setup

- Model: RandomForest(500 trees, class_weight=balanced_subsample) (434 input features). Features: amino-acid & dipeptide composition, residue groups, physico-chemistry.
- Time split: train = reviewed entries first published <= **2018-12-31**; hold-out = published after.
- Filter: `reviewed:true AND fragment:false AND length:[50 TO 2000]`; entries with >1 of the class keywords excluded.
- Classes: kinase (KW-0418), protease (KW-0645), ion_channel (KW-0407), ribosomal_protein (KW-0689), toxin (KW-0800), dna_binding (KW-0238), other (none of the other keywords)
- Note: ion_channel: only 69 holdout entries available (wanted 100).
- Note: Removed 71 duplicate accessions/sequences.

| class             |   train |   hold-out | hold-out recall   |
|-------------------|---------|------------|-------------------|
| dna_binding       |     399 |        100 | 66.0%             |
| ion_channel       |     388 |         68 | 77.9%             |
| kinase            |     377 |        100 | 23.0%             |
| other             |     400 |        100 | 16.0%             |
| protease          |     396 |         99 | 62.6%             |
| ribosomal_protein |     376 |         96 | 84.4%             |
| toxin             |     399 |        100 | 76.0%             |

## Hold-out results (time split)

| method                                                               | accuracy   | balanced acc.   | macro-F1   |
|----------------------------------------------------------------------|------------|-----------------|------------|
| **Random forest on 434 hand-crafted sequence features (this model)** | 56.9%      | 58.0%           | 55.8%      |
| majority class                                                       | 15.1%      | 14.3%           | 3.7%       |
| random (stratified)                                                  | 13.9%      | 14.2%           | 13.9%      |
| logistic regression, composition only (20 AAC)                       | 45.9%      | 46.8%           | 45.7%      |
| 1-NN k-mer homology transfer                                         | 49.2%      | 49.3%           | 49.0%      |

## Leakage check: accuracy vs. similarity to the training set

For each hold-out protein we find its 5 most similar training proteins (3-mer Jaccard), align it to each (local alignment; identical residues / shorter length) and keep the highest identity. A newer entry that is a close homolog of an older one is easy; the low-similarity rows are the honest measure of generalisation to genuinely new proteins.

| 3-mer similarity      |   n | accuracy   |
|-----------------------|-----|------------|
| <0.10 (dissimilar)    | 346 | 55.5%      |
| 0.10-0.20             | 240 | 52.1%      |
| 0.20-0.40             |  47 | 72.3%      |
| >=0.40 (near-homolog) |  30 | 86.7%      |

| identity to nearest train protein   |   n | accuracy   |
|-------------------------------------|-----|------------|
| <30% identity                       | 482 | 50.6%      |
| 30-50%                              |  84 | 66.7%      |
| 50-70%                              |  45 | 73.3%      |
| >=70%                               |  52 | 84.6%      |

- Median nearest-neighbour similarity: 0.097; median identity: 10.7%
- **Low-similarity hold-out (3-mer < 0.10, n=346): accuracy 55.5%, balanced 52.8%**
- **Hold-out with < 30% identity to any training protein (n=482): accuracy 50.6%, balanced 54.2%**

### Redundancy-filtered hold-out

Removing the 139 hold-out proteins that are >= 40% identical to any of their 5 closest training proteins leaves n=524:

| method                 | accuracy   | balanced acc.   | macro-F1   |
|------------------------|------------|-----------------|------------|
| this model             | 52.3%      | 55.7%           | 53.4%      |
| 1-NN homology transfer | 35.9%      |                 |            |

This approximates clustering at 40% identity (as MMseqs2/CD-HIT would) for the train/hold-out boundary. Candidates come from a 3-mer prefilter, so a rare distant homolog could be missed.

## Confusion matrix (rows = true, columns = predicted)

| true \ pred       |   dna_binding |   ion_channel |   kinase |   other |   protease |   ribosomal_protein |   toxin |
|-------------------|---------------|---------------|----------|---------|------------|---------------------|---------|
| dna_binding       |            66 |             3 |        4 |       6 |         10 |                  10 |       1 |
| ion_channel       |             7 |            53 |        0 |       1 |          6 |                   0 |       1 |
| kinase            |            39 |             8 |       23 |      16 |         14 |                   0 |       0 |
| other             |            28 |             7 |       11 |      16 |         32 |                   3 |       3 |
| protease          |            12 |             5 |        8 |       9 |         62 |                   0 |       3 |
| ribosomal_protein |             8 |             0 |        0 |       2 |          1 |                  81 |       4 |
| toxin             |             5 |             7 |        0 |       4 |          3 |                   5 |      76 |

## Top features

| feature            |   importance |
|--------------------|--------------|
| n_tm_segments      |      0.02776 |
| pI                 |      0.02277 |
| aac_C              |      0.01881 |
| charge_per_residue |      0.01865 |
| group_positive     |      0.01805 |
| log10_length       |      0.01603 |
| aromaticity        |      0.01299 |
| gravy              |      0.01207 |
| aac_R              |      0.01068 |
| group_hydrophobic  |      0.00919 |
| aac_F              |      0.00904 |
| dpc_CC             |      0.00885 |
| aac_K              |      0.00825 |
| aac_L              |      0.00796 |
| group_aromatic     |      0.00779 |

## Class patterns

| class             |   n |   median len |   GRAVY |    pI |   frac. TM | enriched residues               |
|-------------------|-----|--------------|---------|-------|------------|---------------------------------|
| dna_binding       | 499 |          389 |  -0.423 |  7.2  |      0.122 | Q (+0.9%), E (+0.9%), R (+0.7%) |
| ion_channel       | 456 |          444 |   0.144 |  7.34 |      0.967 | L (+2.4%), F (+1.9%), I (+1.2%) |
| kinase            | 477 |          476 |  -0.221 |  6.43 |      0.279 | L (+0.8%), E (+0.7%), D (+0.7%) |
| other             | 500 |          391 |  -0.214 |  6.51 |      0.312 | A (+0.6%), D (+0.5%), Q (+0.4%) |
| protease          | 495 |          506 |  -0.281 |  6.58 |      0.535 | S (+0.7%), N (+0.6%), D (+0.5%) |
| ribosomal_protein | 472 |          165 |  -0.491 | 10    |      0.028 | K (+3.5%), R (+3.1%), V (+2.0%) |
| toxin             | 499 |          113 |  -0.208 |  7.34 |      0.707 | C (+5.0%), K (+1.2%), Y (+1.1%) |

## Figures

![confusion_matrix_handcrafted](confusion_matrix_handcrafted.png)
![accuracy_by_similarity_handcrafted](accuracy_by_similarity_handcrafted.png)
![feature_importance_handcrafted](feature_importance_handcrafted.png)
![pca_map_handcrafted](pca_map_handcrafted.png)
![property_boxplots](property_boxplots.png)
