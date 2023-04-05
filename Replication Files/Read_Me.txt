THE SOLVER USED IS KNITROMATLAB. IF THE KNITROMATLAB SOVLER IS UNAVAILABLE, YOU CAN USE FMINCON.
THERE ARE SOME FUNCTIONS USED THAT NEED TO BE ADDED TO ROOT DIRECTORY: kurtosis.m, mvnrnd.m, parfor_progress.m, skewness.m, and matlab2tikz.m

* Growth accounting results of Section 6.3 can be replicated by running

"\Growth Accounting_Klems\Growth_Accounting.m"

* GDP Simulation results of Section 6.1 can be replicated by running (in order)

"GDP Simulatin --- 88 Sector\GDP_Simulation_88sectorKLEMS.m"
"GDP Simulatin --- 88 Sector\GDP_Simulation_88sectorKLEMS_robustness.m"
"GDP Simulatin --- 88 Sector\GDP_Simulation_88sectorKLEMS_robustness_reallocation.m"
"\GDP Simulatin -- 88 Sector\Robustness\Create_Table.m"

The first replicates results with immobile labor and the second replicates results with mobile labor. The last one generates the oil vs. retail picture and well as the oil vs. construction picture. The last makes the robustness tables in Online Appendix 3. 

* GDP Simulation results with adjustment costs of Online Appendix 2 can be replicated by running 

"\Stuck Intermediates and Adjustment Costs\Run_Simulations_with_Adjustment_Costs.m"
"\Stuck Intermediates and Adjustment Costs\Results with Kappa\Create_Table_Adj_costs.m"

The first file runs a given simulation. The second file creates Table 1 in Online Appendix 2. 

