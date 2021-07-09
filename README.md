# MultiplexedNormalization
R code and scripts to run the multiplexed imaging normalization methods, detailed in the "Quantifying and correcting slide-to-slide variation in multiplexed immunofluorescence images" paper. Please see [here](https://colemanrharris.me) to explore the paper.

The `normalization_methods.Rmd` file performs all of the necessary analyses to recreate the analysis performed in the paper - all necessary packages and setup are performed in the markdown file itself (an HTML file generated with `knitr` is included in this repository as well). The analyses performed include:
- Loading the dataset
- Calculating transformations on the desired marker channels
- Implementing the functional data registration algorithm
- Implementing the ComBat algorithm
- Generating Otsu thresholds and metrics as defined in the paper
- Recreating each figure as presented in the paper (incl. random effects modeling and UMAP embedding).

This R markdown file relies on R scripts as included in the repository, that implement the following mechanisms:
- Performing the functional data registration for multiplexed imaging (`all_fda_functions_210709.R`)
- Implementing the ComBat algorithm for multiplexed imaging (`all_combat_functions_210709.R`)
- Calculting Otsu thresolds and metrics (`all_otsu_functions_210709.R`)

Please feel free to reach out to me, coleman.r.harris@vanderbilt.edu, with any questions.
