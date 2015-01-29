
Rationality In Switching Environments (RISE) Toolbox
================
Welcome to RISE!!! For any issue, suggestion or bug
report, please send an email to [junior.maih AT gmail.com](junior.maih@gmail.com)

RISE is an object-oriented Matlab toolbox for solving and estimating nonlinear
dynamic stochastic general equilibirium (henceforward, DSGE) or more generally
Rational Expectations(hereinafter RE), models with switching parameters.

Leading references in the field include various papers by Dan Waggoner and Tao Zha, which
can be found [here](http://www.tzha.net/articles).

In RISE, the switching of the parameters is governed by Markov processes and can be endogenous.

RISE uses perturbation to approximate the nonlinear Markov Switching Rational
Expectations (MSRE) model and solves it using efficient algorithms.

Constant-parameter DSGE/RE models are a special case. RISE also includes the solution
and estimation of
* [Loose commitment](http://journals.cambridge.org/action/displayAbstract?fromPage=online&aid=8686985) (or optimal policy) models
* sticky information models
* hybrid expectations models
* [models in which agents have information about future events](http://www.kansascityfed.org/publicat/events/research/2010CenBankForecasting/Maih_paper.pdf)

Features of the toolbox
-----------------------------------

`*` - Symbolic, automatic and numerical derivatives

`*` - Time Series

`*` - Maximum Likelihood and Bayesian estimation

`*` - Derivative-free optimization algorithms

`*` - Forecasting and conditional forecasting

`*` - DSGE-VAR modeling

`*` - Integrated SVAR modeling with short, long and mixed restrictions (Recently added)

`*` - High Dimensional Model Representation and Monte Carlo Filtering

`*` - Reporting	system

Installation
-----

1. Download the toolbox either as zip file or through github (recommended)
2. add the toolbox path to Matlab and run rise_startup();

Example: Foerster, Rubio-Ramirez, Waggoner and Zha (2013)
---------------------
```RISE
% the example can be found under RISE_toolbox\examples\MarkovSwitching\RamirezWaggonerZha2012
% and the name of the file is frwz_nk.dyn

//Declare the endogenous
endogenous	PAI,Y,R

//Declare the exogenous
exogenous EPS_R

//Declare the parameters: a_tp_1_2 and a_tp_2_1 are automatically recognized by RISE
// as transition probabilities. In this particular case, the name of the markov chain
// is "a"; "tp" stands for transition probability; and in "i_j", "i" denotes the
// the current regime and "j" denotes the regime next period.
parameters a_tp_1_2, a_tp_2_1, betta, eta, kappa, mu,
mu_bar, psi, rhor
sigr

// equations of the model
model
	1-betta*(1-.5*kappa*(PAI-1)^2)*Y*R/((1-.5*kappa*(PAI(+1)-1)^2)*Y(+1)*exp(mu)*PAI(+1));
	
	1-eta+eta*(1-.5*kappa*(PAI-1)^2)*Y+betta*kappa*(1-.5*kappa*(PAI-1)^2)*(PAI(+1)-1)*PAI(+1)/(1-.5*kappa*(PAI(+1)-1)^2)
	-kappa*(PAI-1)*PAI;

	(R(-1)/steady_state(R))^rhor*(PAI/steady_state(PAI))^((1-rhor)*psi)*exp(sigr*EPS_R)-R/steady_state(R);

	
// the steady state
steady_state_model(unique,imposed)
    PAI=1;
    Y=(eta-1)/eta;
    R=exp(mu_bar)/betta*PAI;

	
// parameterization
parameterization
	a_tp_1_2,1-.9; 
	a_tp_2_1,1-.9;
	betta, .99;
	kappa, 161;
	eta, 10;
	rhor, .8;
	sigr, 0.0025;
	mu_bar,0.02; 
	mu(a,1), 0.03; // value assumed by parameter "mu" (controled by chain a) in regime 1
	mu(a,2), 0.01;
	psi(a,1), 3.1;
	psi(a,2), 0.9; // value assumed by parameter "psi" (controled by chain a) in regime 2

```

Running the example
---------------------
```matlab

%% housekeeping
clear all;
close all;
clc;

%%  load the model and its parameterization

frwz=rise('frwz_nk');

%% Solving the model

frwz=solve(frwz);

%% print the solution

print_solution(frwz)

```
Documentation
---------------------
* This publicly available version of RISE is largely work in progress and there
is a documentation (PDF, HTML, etc) but that too is work in progress.
* RISE is very intuitive as you will quickly realize when you run the examples in RISE_toolbox\examples\
* Should you want to report a bug, make a suggestion or need help, please send an email to [junior.maih AT gmail.com](junior.maih@gmail.com) 