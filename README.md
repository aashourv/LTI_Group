# LTI_Group

LTI_GROUP estimates the parameters of a first-order multivariate AR model 
for each subject's multivariate time series.

Syntax:
  [A_opt, J_opt] = LTI_Group(V)

Input:
  V - A cell array where each cell V{subject_i} is an (n by t) matrix 
  containing the time series data for subject i. Rows of V{subject_i} 
  represent variables, and columns represent observations.

Output:
  A_opt - Estimated coefficient matrix (n by n).
  J_opt - The value of the cost function at the optimal A.

Description:
  This function estimates the parameters of a first-order multivariate AR 
  model for each subject's multivariate time series. The model is defined 
  as:

      x(:,k) = A*x(:,k-1) + noise,

  where x is the state vector and A is the system matrix to be estimated. 
  The function provides least squares estimates of the coefficient matrix 
  A (n*n) and assumes that the AR process has zero mean.

  The optimization problem is defined as minimizing the sum of squared 
  errors across all subjects:

      J(A) = \sum_{i=1}^{N} \sum_{k=1}^{T} (x_i[k] - \hat{x}_i[k])^2,

  where x_i[k] is the observed data for subject i at time k, and 
  \hat{x}_i[k] is the model prediction based on the LTI system with 
  system matrix A. The optimization is performed using iterative 
  algorithms such as the quasi-Newton method.

  To find the optimal matrix A, we used the fminunc function as follows:

      [A_opt, J_opt] = fminunc(@(A) cost_function(A, V), A0, options);

  - fminunc: This is a function used for unconstrained optimization. It 
    finds the minimum of an unconstrained multivariable function.
  - @(A) cost_function(A, V): This creates an anonymous function that 
    takes A as an input and calls the cost_function with A and V. The 
    cost_function calculates the sum of squared errors between the 
    observed data and the model predictions.
  - A0: This is the initial guess for the system matrix A. It is an 
    identity matrix of appropriate size.
  - options: These are optimization options specified by the optimoptions 
    function. In this case, 'Display' is set to 'iter-detailed' to provide 
    detailed information at each iteration of the optimization process.
  - A_opt: This is the output of fminunc, representing the optimal system 
    matrix A that minimizes the cost function.
  - J_opt: This is the value of the cost function at the optimal A.

Example:
  V = {randn(5, 100), randn(5, 100)};
  [A_opt, J_opt] = LTI_Group(V);

See also:
  FMINUNC, OPTIMOPTIONS

Please contact Arian Ashourvan (ashourvan@ku.edu) with any questions regarding this code.
  
Authors:
  Arian Ashourvan, University of Kansas
  Sérgio Pequito, Universidade de Lisboa

Date:
  July 2024
