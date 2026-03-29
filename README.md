# RISE Toolbox

**Rationality In Switching Environments (RISE)** is an object-oriented toolbox written in MATLAB, primarily designed for the modelling of regime-switching dynamic stochastic general equilibrium (DSGE) models. The transition probabilities governing the switching mechanism can be exogenous or endogenous. Constant-parameter models, which are special cases, are also fully supported.

> **Legacy release:** The previous public version of RISE (circa 2019) is preserved as a tagged release [`v2019`](../../releases/tag/v2019) for users who need backward compatibility. Note that the legacy version is no longer supported.

---

## Installation

1. Clone or download this repository.
2. Add the toolbox path to MATLAB and run `rise_startup()`.

---

## Features

### DSGE Modelling

- Perturbation solution up to 5th order for regime-switching DSGE models
- Endogenous and exogenous switching mechanisms
- Loose commitment policies with and without occasionally binding constraints — see [Debortoli, Maih and Nunes (2012)](http://journals.cambridge.org/action/displayAbstract?fromPage=online&aid=8686985)
- 5th-order perturbation solution of regime-switching optimal-policy DSGE models: commitment, discretion, and stochastic replanning, with and without occasionally binding constraints
- Deterministic solution of nonlinear regime-switching DSGE models
- Conditional forecasting in linear and nonlinear models
- Models in which agents have information about future events
- Automatic translation of files from other toolboxes
- Powerful macro-language

### VAR Modelling

- Reduced-form VARs with and without regime switching
- Panel reduced-form VARs with and without regime switching
- Structural VARs with and without regime switching
- Proxy structural VARs with and without regime switching

### DSGE-VAR

- Estimation and simulation of DSGE-VAR models

### Estimation

- Maximum likelihood estimation
- Bayesian estimation
- Indirect inference

### Uncertainty Quantification

- Global sensitivity analysis
- High-dimensional model representation

### Global Optimisation

- Derivative-free stochastic global optimisation algorithms

### Other

- Time series and data management
- Reporting system

---

## Key References

The main methodological paper describing the regime-switching problem and its solution as implemented in RISE:

- [Junior Maih (2015): "Efficient perturbation methods for solving regime-switching DSGE models"](http://www.norges-bank.no/en/Published/Papers/Working-Papers/2015/12015/)

Other useful references include various papers by Dan Waggoner and Tao Zha, available [here](http://www.tzha.net/articles).

---

## Selected Applications

- [Hashimzade, Kirsanov, Kirsanova and Maih (2026): "Filtering and Smoothing in State-Space Models with Multiple Regimes" forthcoming in the *Journal of Business and Economic Statistics*]
- [Yoosoon Chang, Junior Maih and Fei Tan (2021): "Origins of monetary policy shifts: A new approach to regime switching in DSGE models", Journal of Economic Dynamics and Control
Volume 133, December 2021, 104235](https://www.sciencedirect.com/science/article/abs/pii/S0165188921001706)
<!-- - [Ryo Kato, Junior Maih and Shin-Ichi Nishiyama (2022): "Trend inflation in Japanese pre-2000s: A Markov-switching DSGE estimation"](http://www.econ.kobe-u.ac.jp/RePEc/koe/wpaper/2022/2212.pdf) -->
- [Ragna Alstadheim, Hilde C. Bjørnland and Junior Maih (2021): "Do central banks respond to exchange rate movements? A Markov-switching structural investigation of commodity exporters and importers", Energy Economics
Volume 96, April 2021, 105138](https://www.sciencedirect.com/science/article/pii/S0140988321000438?via%3Dihub)
- [Junior Maih, Falk Mazelis, Roberto Motto and Annukka Ristiniemi (2021): "Asymmetric monetary policy rules for the euro area and the US", Journal of Macroeconomics
Volume 70, December 2021, 103376](https://www.sciencedirect.com/science/article/abs/pii/S0164070421000756)
- [Andrew Binning, Hilde C. Bjørnland and Junior Maih (2019): "Is monetary policy always effective? Incomplete interest rate pass-through in a DSGE model"](https://www.norges-bank.no/aktuelt/nyheter-og-hendelser/Signerte-publikasjoner/Working-Papers/2019/222019/)
- [Hilde C. Bjørnland, Vegard H. Larsen and Junior Maih (2018): "Oil and macroeconomic (in)stability", American Economic Journal: Macroeconomics
vol. 10, no. 4, October 2018
(pp. 128–51)](https://www.aeaweb.org/articles?id=10.1257/mac.20150171)
- [Andrew Binning and Junior Maih (2017): "Modelling occasionally binding constraints using regime-switching"](https://www.norges-bank.no/en/news-events/news-publications/Papers/Working-Papers/2017/232017/)
- [Andrew Binning and Junior Maih (2016): "Forecast uncertainty in the neighborhood of the effective lower bound: How much asymmetry should we expect?"](https://www.norges-bank.no/en/news-events/news-publications/Papers/Working-Papers/2016/132016/)
- [Farooq Akram, Andrew Binning and Junior Maih (2016): "Joint prediction bands for macroeconomic risk management"](https://www.norges-bank.no/en/news-events/news-publications/Papers/Working-Papers/2016/72016/)
- [Andrew Binning and Junior Maih (2016): "Implementing the zero lower bound in an estimated regime-switching DSGE model"](https://www.norges-bank.no/en/news-events/news-publications/Papers/Working-Papers/2016/32016/)
- [Andrew Binning and Junior Maih (2015): "Applying flexible parameter restrictions in Markov-switching vector autoregression models"](https://www.norges-bank.no/en/news-events/news-publications/Papers/Working-Papers/2015/172015/)
- [Andrew Binning and Junior Maih (2015): "Sigma point filters for dynamic nonlinear regime switching models"](https://www.norges-bank.no/en/news-events/news-publications/Papers/Working-Papers/2015/102015/)

---

## Extended Version

An extended version of RISE is available to contributors and participants in RISE events. Contact the author for details.

## Getting Help and Reporting Issues

For bug reports, suggestions, or questions, please open an issue on the [GitHub repository](https://github.com/jmaih/RISE_toolbox/issues) or send an email to [junior.maih AT gmail.com](mailto:junior.maih@gmail.com).