# Ecological Behaviour of Minimal Reactomes (Modified Pipeline)

This repository contains the upgraded MATLAB pipeline for extracting Minimal Reactomes from 52 AGORA2 gut microbiome models and analyzing their pairwise ecological interactions under two dietary conditions.

## Pipeline Overview

### 1. Diet Application (`applyDiet.m`)
- Applies two physiological diets: **High Fiber** and **Western** (from VMH database).
- Standardises metabolite nomenclature (`EX_h2`, `EX_glc_D`).

### 2. Minimal Reactome Extraction (`Minimal_Reactomes.m`)
- **Wild-Type Rescue:** Models that cannot grow on a diet are rescued using `minimalMedium()` from the COBRA Toolbox. The lacking media components are supplemented and the model is flagged as `Supplemented`.
- **Three Reduction Scenarios:**
  - **Scenario 1 (`cut1`):** `GrowthRateCutoff = 1.0` — near-WT growth required (strict).
  - **Scenario 2 (`cut05`):** `GrowthRateCutoff = 0.5` — 50% WT growth sufficient (relaxed, larger genome retained).
  - **Scenario 3 (`ext`):** `GrowthRateCutoff = 0.5` + all exchange and transport reactions preserved in the elite list (ensures community interaction capability).

### 3. Pairwise Community Analysis (`Pairwise_Community.m`)
- Communities are grown on the **original, unsupplemented diets** (not the supplemented diet used during rescue).
- Evaluates WT and all three MRM variants in pairwise combinations (1,326 pairs).
- Computes FVA-based growth ranges at the community optimum for robust interaction classification.

### 4. Visualization (`Plot1.m` – `Plot7.m`)
- **Plot 1:** Interaction types stacked bar chart (WT + 3 MRM cuts × 2 diets)
- **Plot 2:** Monoculture vs community growth scatter (4×2 grid)
- **Plot 3:** Interaction composition as normalised percentages
- **Plot 4:** 52×52 heatmap (upper triangle = WT, lower triangle = MRM)
- **Plot 5:** Sankey transition diagram (WT → MRM interaction shifts)
- **Plot 6:** Phylum-pair interaction composition (scaled pie charts)
- **Plot 7:** Jaccard similarity vs interaction type (tertile bins)
- Helper functions: `calcAlpha.m`, `classifyInteraction.m`

## Hypothesis
- If genome is reduced, interaction is more (organisms depend on partners).
- If genome is big, organisms are self-sustaining, so interactions are less.
- If genome is very small, the potential for interaction (exchanges) might also get reduced.
- Changing the growth cutoff changes the genome size, which changes the interaction type.

## Execution
Run the scripts in order using MATLAB R2019b+ with the COBRA Toolbox and Gurobi:
1. `Minimal_Reactomes.m`
2. `Pairwise_Community.m`
3. `Plot1.m` through `Plot7.m` (require output CSVs from steps 1–2)

## Requirements
- MATLAB R2019b or higher
- Gurobi Optimizer
- COBRA Toolbox v3.0+

## Author
**Avinash Chauhan**
GitHub: [@avinash-chauhan-oss](https://github.com/avinash-chauhan-oss)
