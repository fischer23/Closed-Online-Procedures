# Closed-Online-Procedures

This repository contains all the code used for the simulations in the paper "The Online Closure Principle" (https://doi.org/10.48550/arXiv.2211.11400). The procedures (ADDIS-Spending and closed ADDIS-Spending) can be found in procedures.R and the code that generates the plots and tables in the Plot_XXX.R files.  The IMPC_data.csv file contains 5000 observations from a full dataset available on a Zenodo repository organized by Robertson et al. (2018) (https://doi.org/10.5281/zenodo.2396572).

In the aforementioned paper, we introduced a new online closure principle including a condition under which the resulting closed testing procedure is indeed an online multiple testing procedure. We showed that each online procedure with FWER control can be derived by this online closure principle and demonstrated how short-cuts of these online closed procedures can be obtained. This was used to define new online multiple testing procedures with FWER control. In particular, we proposed the closed ADDIS-Spending as an uniform improvement of the ADDIS-Spending procedure (Tian and Ramdas, 2021). 


Robertson, D. S., Wildenhain, J., Javanmard, A., & Karp, N. A. (2018). Supporting data for "onlineFDR: an R package to control the false discovery rate for growing data repositories" (1.2) [Data set]. Zenodo. https://doi.org/10.5281/zenodo.2396572
