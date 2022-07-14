# Predicting the Heat of Formation for Materials with 1–3 Elements
This was the capstone project for my Master's degree in Data Analytics. It covers a basic application of machine learning in materials informatics (MI): the use of existing materials properties to predict another property for use when, for example, assessing new materials for research and development.

## Data
The raw data were downloaded from the [Computational Materials Repository (CMR)](https://cmr.fysik.dtu.dk/oqmd123/oqmd123.html). These data include the chemical formula, the volume of the unit cell (in cubic angstroms), the sum of the atomic masses in the unit cell (in grams per mole), the number of atoms in the chemical formula, and the heat of formation (in electron-volts per atom).

Additional features related to the elements in each compound were created using [Magpie](https://bitbucket.org/wolverton/magpie/).

## Models
For elemental compounds, the heat of formation is zero according to thermodynamics. Thus, the overall model will predict that any new one-component material will have a heat of formation of zero. For two- and three- component materials, multiple different regression and tree-based models were built to predict the heat of formation. The best model uses boosted decision trees and achieves a mean absolute error of 68 meV/atom, exceeding a literature benchmark of 96 meV/atom derived from density functional theory calculations and experiments. The model performance also compares favorably to other models in the literature (for details, see the final report)

## Citation
If you find this work useful, please cite it as

I. P. Steinke, "Predicting the Heat of Formation of One-, Two-, and Three-Component Materials Using Machine Learning," M.S. Capstone Project - Final Report, Slippery Rock University, Slippery Rock, PA, 2021.
