

# TACO: Transitional Analysis through Clymo Optimization

The TACO (Transitional Analysis through Clymo Optimization) ers a novel Bayesian approach to understanding the complex carbon dynamics within peatland ecosystems, extending beyond traditional models to provide a more nuanced view of carbon flux and decay processes. Inspired by the classic Clymo model, this package introduces a sophisticated Bayesian framework capable of accommodating variable carbon influxes even in the deepest sediment layers and allowing for smooth transitions between different system states within peatlands.


### Installation

Download the Taco.R and the twalk.R files. soruce the file:

```R
source(\PATH\TO\Taco.R)

```


### Features

- **Bayesian Inference:** Utilize Bayesian statistics to make accurate inferences about peatland carbon dynamics.
- **Weighted Function:** A unique weighted function facilitates the modeling of transitions between different peatland systems, reflecting changes in carbon decomposition rates.
- **Flexible Modeling:** Accommodates variable carbon influxes, overcoming limitations of traditional models that fail to capture these variations, especially in deeper sediment layers.
- **Peatland Management Implications:** Offers insights into peatland management and the global carbon cycle, supporting research and conservation efforts.

### Usage

The package includes functions for constructing the Bayesian model, performing inference, and analyzing peatland carbon dynamics. Example usage:

```R
setwd('~/Documents/Taco/')
source('TACO.R')
tac <- Taco(CORE_NAME , '~/Documents/Taco/')
```

### Documentation

_Working on it_


### License

This package is released under the MIT License. See the LICENSE file for more details.

### Citation

Manuscript in process

### Contact

For any questions or feedback, please contact Marco A Aquino-Lopez at aquino@cimat.mx.