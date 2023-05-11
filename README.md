# Bayesian Weighted Clymo

This is a variation of the acrotelm and catotelm modelling approach. In the traditional approach a "hard" boundary is chosen and two different Clymo approaches are fitted to the data. In order to allow for a more realistic transition, this approach considers a weigh function ($w(x)$), which range is between 0 and 1 and domain is within the length of the core. In this model, two Clymo models are assumed to affect the core and then weighted by $w(x)$. Function $w(x)$ is assumed to be monotonic in order to coincide with the idea that for peat which is below the water-table only ca totem effects take effect. 

## Results

The outputs of this method are the posterior distributions of the parameters of both Clymo models. 

![Posterior Distributios of both Climo models](./Figures/Posterior_parameters.pdf)

Because the model also has to infer the $w(x)$ function the posterior distribution of the parameters of the function can also be obtained.


![Posterior Distributios of the $w(x)$ function](./Figures/Bon_posterior.pdf)

Lastly, this posterior parameters can be use to plot the posterior transition function. 

![Posterior transition function](./Figures/Cat_limit.pdf)




