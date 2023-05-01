# Simulation code for transporting subgroup treatment effects

`Reference`: Paper title and link

### `Data_Generation` 
contains parameter specifications and functions to generate data and the true values.

The outcome (Y) is continuous, generated by non-linear functions ofcovariates (X), treatment (A) and their interactions. We propose to use Generalized additive model for estimating the outcome model. 

Other models can be correctly estimated by multi normial or generalized linear regression.

### `Internal_External_Estimates_DR` 
contains the code to show 
1. biases, 
2. coverage of simultaneous confidence bands, 
3. coverage of pointwise confidence intervals,
4. theoretical standard deviation (based on the influence function) and 
5. Monte-Carlo standard deviation 

of each of the subgroup effects using the propsed doubly robust estimator for estimating $\phi_{a=1,s=1}(\tilde{X}=\tilde{x})$ and $\psi_{a=1}(\tilde{X}=\tilde{x}), \forall \tilde{x}=1,\dots,5$ with the specified $n$ (sum of sample sizes of internal and extenal population) and $n_m$ (sample size of internal population). 

The results are show in plots. They show the consistency of the proposed doubly robust estimator, and the relationships between confidence intervla and simutaneous confidence bands.

### `Internal_Mis` 
contains the code to show the biases and standard deviations of 
1. (Doubly robust estimator) $\widehat \phi_{a=1,s=1}(\tilde{X}=3)$, 
2. (G computation-based estimator) $\widehat \phi^{g}_{a=1,s=1}(\tilde{X}=3)$, and 
3. (IPTW-based estimator) $\widehat \phi^{w}_{a=1,s=1}(\tilde{X}=3)$

when

1. Both outcome and propensity score and source models are correctly specified,
2. Only the outcome model is correctly specified,
3. Only the propensity score and source models are correctly specified,
4. None of the models are correctly specified.

with the specified $n$ (sum of sample sizes of internal and extenal population) and $n_m$ (sample size of internal population).  

The results are show in plots. They show the doubly robustness property of teh estimator.

### `External_Mis` 
contains the code to show the biases and standard deviations of 
1. (Doubly robust estimator) $\widehat \psi_{a=1}(\tilde{X}=3)$, 
2. (G computation-based estimator) $\widehat \psi^{g}_{a=1}(\tilde{X}=3)$, and 
3. (IPTW-based estimator) $\widehat \psi^{w}_{a=1}(\tilde{X}=3)$

when

1. Both outcome and propensity score and source models are correctly specified,
2. Only the outcome model is correctly specified,
3. Only the propensity score and source models are correctly specified,
4. None of the models are correctly specified.

with the specified $n$ (sum of sample sizes of internal and extenal population) and $n_m$ (sample size of internal population). 

The results are show in plots.  They show the doubly robustness property of teh estimator.

### `Rates` 
contains the code for estimating the RMSE with the varying convergence rates of the
1. (Doubly robust estimator) $\widehat \phi_{a=1,s=1}(\tilde{X}=1)$, 
2. (G computation-based estimator) $\widehat \phi^{g}_{a=1,s=1}(\tilde{X}=1)$

with the specified $n_m$ (sample size of internal population). 

The results are show in plots. They show the rate robustness of the estimator.


