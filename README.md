# Replication Files for Universal Difference-in-Differences for Causal Inference in Epidemiology (Tchetgen Tchetgen, Park, Richardson, _Epidemiology_, 2023+) 

This Github repository contains replication files for [Universal Difference-in-Differences (UDiD) for Causal Inference in Epidemiology (Tchetgen Tchetgen, Park, Richardson, 2023)](https://arxiv.org/abs/2302.00840 "UDiD").


## Data

The dataset contains birth rate information from 603 municipalities in two states of Brazil, Pernambuco and Rio Grande do Sul. 
Municipality-level birth rates were measured in 2014 and 2016, before and after the 2015 Zika virus outbreak.
More details on the source of the dataset are given below:
* Pre- and Post-treatment Outcomes, Treatment, and log population: zika_Table2.tab in  https://dataverse.harvard.edu/dataset.xhtml?persistentId=doi:10.7910/DVN/ENG0IY. 
See [Taddeo, Amorim, Aquino (2022)](https://www.intlpress.com/site/pub/pages/journals/items/sii/content/vols/0015/0004/a001/index.php?mode=ns "ZB") for details.
* log population density and proportion of female: https://www.ibge.gov.br/en/statistics/social/income-expenditure-and-consumption/18391-2010-population-census.html?=&t=resultados


## Code

* Brazil_Zika.R replicates the main analysis in Section 4 and Assessment of Covariate Distribution Invariance in Section A.10 of the Appendix.
* UDID.R contains functions used in Brazil_Zika.R

## References

Taddeo, Amorim, Aquino (2022) **Causal Measures Using Generalized Difference-in-difference Approach with Nonlinear Models**, _Statistics and Its Interface_, 15(4):399-413 [[link](https://www.intlpress.com/site/pub/pages/journals/items/sii/content/vols/0015/0004/a001/index.php?mode=ns "ZB")]

Tchetgen Tchetgen, Park, Richardson (2023+) **Universal Difference-in-Differences for Causal Inference in Epidemiology**, _Epidemiology (In Press)_ [[arXiv link](https://arxiv.org/abs/2302.00840 "UDiD")]
