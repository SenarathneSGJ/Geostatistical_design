# Geostatistical_design

1. Example 1:
1.1 Obtaining Optimal Designs
 The r-code to obtain designs using three loss functions are available in the following folder: Example1/Rcodes/design_selection.
 It is required to install the r-package MaternEx1_1.0.tar.gz to run this example.
 The designs for this example can be obtained by running the main r-script called main_Clayton.R inside the folder.
 Inside the file main_Clayton.R, there are different options for users to run this code. They are:
  o d_no = number of design points in the design
  o utility = negative value of the loss function
  o Dependence = spatial dependency structure of the two responses
We have run this code with two different values for d_no (5 and 10) and three different values for Dependence (0.2, 0.5 and 0.8) based on the three loss functions. The resulting designs are saved in the following folder: Example1/Selected_designs.

1.2 Design Evaluation
 After designs have been determined for each loss function, they can be evaluated based on each design objective.
 For this purpose, the Example1_design_evaluation.R function can be used. Here, two parameters, d_no and Dependence should be set before running this code.
 Design evaluation results are saved in the folder Example1/Simulation_results/ Design_evaluation.

1.3 Compare two approximations
 To compare the two approximations used in the paper for evaluating the prediction loss values, the function main.R can be used. This file is location in folder: Example1/Rcodes/ Comparing_two_approximations.

1.4 Simulation Results
 To obtain the plots and tables shown for Example 1 we used the results obtained in Sections 1.2 and 1.3 above.
 The results obtained in Sections 1.2 and 1.3 are available in Example1/Simulation _results.
 Example1_plots.R in the Master folder can be used to plot the results for Example 1.

2. Example 2:
2.1 Obtaining Optimal Designs
 The r-code to obtain designs using three loss functions are available in the following folder: Example2/Rcodes/design_selection.
 It is required to install the r-package AirQualityRcpp_1.0.tar.gz to run this example.
 The designs for this example can be obtained by running the main script called main_Clayton.R inside the folder.
 Inside the main_Clayton.R, there are two options for users to run this code. They are,
  o d_no = number of design points in the design
  o utility = negative value of the loss function
 We have run this code with four different values for d_no (5, 7, 10 and 15) based on the three loss functions. The resulting designs are saved in the following folder: Example1/Selected_designs.

2.2 Design Evaluation
 After designs have been determined for each loss function, they can be evaluated based on each design objective.
 For this purpose, the function called Example2_design_evaluation.R can be used. Here, two parameters, d_no and Dependence should be set before running this code.
 Design evaluation results are saved in the folder: Example2/Simulation _results/ Design_evaluation.

2.3 Compare two approximations
To compare the two approximations used in the paper for evaluating the prediction loss values, the main.R function in the following folder can be used: Example2/Rcode/Comparing_two_approximations.

2.4 Simulation results
 To obtain the plots and tables shown for Example 2 we used the results obtained in Sections 2.2 and 2.3 above.
 The intermediate results obtained in Sections 2.2 and 2.3 are saved in Example2/Simulation_results.
 We use the function Example2_plots.R available in the Master folder to obtain the results related to Example 2.

3. Design Comparison (Supplementary material):
 To obtain designs from cluster 3 in Example2, we used the same set of code discussed in Section 2.1 above.
 Here, the function main_Clayton_C3.R is used for design selection.
 To compare the efficiencies of the two approximations to the prediction loss function, the same set of r-code discussed in Sections 1.3 and 2.3 can used. Here, the main_efficiency.R functions available in the respective folders can be used.
 The figures and tables available in the Supplementary material can be obtained using the function supplementary_plots.R available in the Master folder.
