# Closed-Online-Procedures

This repository contains all the code used for the simulations in the paper "The Online Closure Principle". The procedures (ADDIS-Spending and closed ADDIS-Spending) can be found in OCP_Procedures.R. The file plot_creator.R applies plot_generating_function.R appropriately to generate the results for the different simulation settings considered in the paper. The real data applications can be reproduced by plot_IMPC.R and plot_platform_RECOVERY.R. The IMPC_data.csv file contains 5000 observations from a full dataset available at a Zenodo repository organized by Robertson et al. (2018) (https://doi.org/10.5281/zenodo.2396572).

In the aforementioned paper, we introduced a new online closure principle including a condition under which the resulting closed testing procedure is indeed an online multiple testing procedure. We showed that each online procedure with FWER control can be derived by this online closure principle and demonstrated how short-cuts of these online closed procedures can be obtained. This was used to define new online multiple testing procedures with FWER control. In particular, we proposed the closed ADDIS-Spending as an uniform improvement of the ADDIS-Spending procedure (Tian and Ramdas, 2021). 

### References

Fischer, L., Roig, M. B. and Brannath, W. (2024). The online closure principle. The Annals of Statistics.

Robertson, D. S., Wildenhain, J., Javanmard, A. and Karp, N. A. (2019). onlineFDR: An R package to control the false discovery rate for growing data repositories. Bioinformatics, 35, 4196–4199.

Tian, J. and Ramdas, A. (2021). Online control of the familywise error rate. Statistical Methods in Medical Research, 30, 976–993.
