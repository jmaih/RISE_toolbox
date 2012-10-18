%% housekeeping
clear all
clc
close all
%% start rise
setpaths
%% "rise" the model
% dynare uses different parameters for baseline calibration and for
% estimation. In rise, it is not allowed to to that. The baseliine
% calibration will be used as starting values for estimation.
% However, thanks to the existence of two separate parameter objects, one
% can change the parameterization in one of the objects.
% THE BOTTOM LINE IS: 
% 1- in the parameterization block, the first column corresponds to both
% the baseline calibration and the starting values for estimation.
% 2- the rise_param object inherits those start values. But now we can
% change them to the original dynare baseline calibration (after reading
% the model)
% that. And so we have to have 2 separate  
% model files in order to do both calibration and estimation

fs2000=rise('fs2000_rise');
%% pushing the particular vector (different from the initial estimation point)
% Note we have to declare the regime of the parameter

calibration={'name', 'value','regime'
	'alp'  , 0.330 , 1
	'bet'  , 0.990 , 1
	'gam'  , 0.003 , 1
	'mst'  , 1.011 , 1
	'rho'  , 0.700 , 1
	'psi'  , 0.787 , 1
	'del'  , 0.020 , 1
	'sig_a', 0.014 , 1
	'sig_m', 0.005 , 1
                    };
fs2000=set_parameters(fs2000,calibration); 

%% solve the model and 
fs2000=solve(fs2000);

%% print the solution 

fs2000.print_solution

disp('This should be identical to the dynare solution (see folder dynare_version)')

%% lets do some estimation: we need the data
% get the data from the dynare_version folder

run dynare_version\fsdat_simul

% rise needs data to be passed as a rise_time_series object with dates and
% names of the variables. So we proceed to constructing the database

% we use the same start date as in the Schorfheide paper even though the
% data are not the same
startdate='1950q1';
database=rise_time_series(startdate,[gp_obs,gy_obs],{'gp_obs','gy_obs'});

%% pass the data to the rise object

% 192 observations are used in dynare
test=rise_date(startdate);
end_date= test.observation_2_date(192);
% the loglinear option of dynare implies:
% 1- we have to take the log of our data. 
% 2- we have to exponentiate the corresponding variables in the rise model
% file. Dynare does this by taking the log of the steady state during
% estimation and this does not ring too transparent to me.

% In rise, we can just take the log of the whole database. This is  
% what we do in passing the data. But we also need to tell rise_time_series
% that the logged variables should have the same names as the original
% ones. In order to do that, we add set the flag to true when taking the
% log.
fs2000=set_options(fs2000,'data',log(database,true),'estim_end_date',end_date);

%% estimate the model

% fs2000=estimate(fs2000);%,'debug',true

profile on
fs2000=estimate(fs2000);%,'optimizer',@csminwellwrap,'debug',true
profile off
profile viewer