%
% Status : main Dynare file 
%
% Warning : this file is generated automatically by Dynare
%           from model file (.mod)

clear all
clear global
tic;
global M_ oo_ options_ ys0_ ex0_ estimation_info
options_ = [];
M_.fname = 'constant';
%
% Some global variables initialization
%
global_initialization;
diary off;
logname_ = 'constant.log';
if exist(logname_, 'file')
    delete(logname_)
end
diary(logname_)
M_.exo_names = 'EZ';
M_.exo_names_tex = 'TFP shock';
M_.exo_names = char(M_.exo_names, 'EBETA');
M_.exo_names_tex = char(M_.exo_names_tex, 'Discount shock');
M_.endo_names = 'W';
M_.endo_names_tex = 'Real wages';
M_.endo_names = char(M_.endo_names, 'N');
M_.endo_names_tex = char(M_.endo_names_tex, 'Labor Hours');
M_.endo_names = char(M_.endo_names, 'C');
M_.endo_names_tex = char(M_.endo_names_tex, 'Consumption');
M_.endo_names = char(M_.endo_names, 'R');
M_.endo_names_tex = char(M_.endo_names_tex, 'Nominal interest rate');
M_.endo_names = char(M_.endo_names, 'BETA');
M_.endo_names_tex = char(M_.endo_names_tex, 'Discount factor');
M_.endo_names = char(M_.endo_names, 'PAI');
M_.endo_names_tex = char(M_.endo_names_tex, 'Inflation rate');
M_.endo_names = char(M_.endo_names, 'Z');
M_.endo_names_tex = char(M_.endo_names_tex, 'TFP process');
M_.endo_names = char(M_.endo_names, 'Y');
M_.endo_names_tex = char(M_.endo_names_tex, 'Output');
M_.endo_names = char(M_.endo_names, 'PSI');
M_.endo_names_tex = char(M_.endo_names_tex, 'Real Marginal Cost');
M_.endo_names = char(M_.endo_names, 'RK');
M_.endo_names_tex = char(M_.endo_names_tex, 'Rental Rate');
M_.endo_names = char(M_.endo_names, 'I');
M_.endo_names_tex = char(M_.endo_names_tex, 'Investment');
M_.endo_names = char(M_.endo_names, 'K');
M_.endo_names_tex = char(M_.endo_names_tex, 'Capital');
M_.endo_names = char(M_.endo_names, 'RTAYLOR');
M_.endo_names_tex = char(M_.endo_names_tex, 'RTAYLOR');
M_.param_names = 'chi';
M_.param_names_tex = 'chi';
M_.param_names = char(M_.param_names, 'eta');
M_.param_names_tex = char(M_.param_names_tex, 'eta');
M_.param_names = char(M_.param_names, 'sigma_');
M_.param_names_tex = char(M_.param_names_tex, 'sigma\_');
M_.param_names = char(M_.param_names, 'betass');
M_.param_names_tex = char(M_.param_names_tex, 'betass');
M_.param_names = char(M_.param_names, 'delta');
M_.param_names_tex = char(M_.param_names_tex, 'delta');
M_.param_names = char(M_.param_names, 'alpha');
M_.param_names_tex = char(M_.param_names_tex, 'alpha');
M_.param_names = char(M_.param_names, 'varphi');
M_.param_names_tex = char(M_.param_names_tex, 'varphi');
M_.param_names = char(M_.param_names, 'theta');
M_.param_names_tex = char(M_.param_names_tex, 'theta');
M_.param_names = char(M_.param_names, 'paistar');
M_.param_names_tex = char(M_.param_names_tex, 'paistar');
M_.param_names = char(M_.param_names, 'phi_pai');
M_.param_names_tex = char(M_.param_names_tex, 'phi\_pai');
M_.param_names = char(M_.param_names, 'phi_y');
M_.param_names_tex = char(M_.param_names_tex, 'phi\_y');
M_.param_names = char(M_.param_names, 'zbar');
M_.param_names_tex = char(M_.param_names_tex, 'zbar');
M_.param_names = char(M_.param_names, 'rhoz');
M_.param_names_tex = char(M_.param_names_tex, 'rhoz');
M_.param_names = char(M_.param_names, 'rhobeta');
M_.param_names_tex = char(M_.param_names_tex, 'rhobeta');
M_.param_names = char(M_.param_names, 'sigmabeta');
M_.param_names_tex = char(M_.param_names_tex, 'sigmabeta');
M_.param_names = char(M_.param_names, 'sigmaz');
M_.param_names_tex = char(M_.param_names_tex, 'sigmaz');
M_.param_names = char(M_.param_names, 'nu');
M_.param_names_tex = char(M_.param_names_tex, 'nu');
M_.param_names = char(M_.param_names, 'g_over_y');
M_.param_names_tex = char(M_.param_names_tex, 'g\_over\_y');
M_.param_names = char(M_.param_names, 'rbar');
M_.param_names_tex = char(M_.param_names_tex, 'rbar');
M_.param_names = char(M_.param_names, 'k_over_y');
M_.param_names_tex = char(M_.param_names_tex, 'k\_over\_y');
M_.param_names = char(M_.param_names, 'n_over_y');
M_.param_names_tex = char(M_.param_names_tex, 'n\_over\_y');
M_.param_names = char(M_.param_names, 'i_over_y');
M_.param_names_tex = char(M_.param_names_tex, 'i\_over\_y');
M_.param_names = char(M_.param_names, 'c_over_y');
M_.param_names_tex = char(M_.param_names_tex, 'c\_over\_y');
M_.exo_det_nbr = 0;
M_.exo_nbr = 2;
M_.endo_nbr = 13;
M_.param_nbr = 23;
M_.orig_endo_nbr = 13;
M_.aux_vars = [];
M_.Sigma_e = zeros(2, 2);
M_.H = 0;
options_.block=0;
options_.bytecode=0;
options_.use_dll=0;
erase_compiled_function('constant_static');
erase_compiled_function('constant_dynamic');
M_.lead_lag_incidence = [
 0 4 0;
 0 5 0;
 0 6 17;
 0 7 0;
 1 8 18;
 0 9 19;
 2 10 0;
 0 11 20;
 0 12 0;
 0 13 21;
 0 14 0;
 3 15 0;
 0 16 0;]';
M_.nstatic = 6;
M_.nfwrd   = 4;
M_.npred   = 2;
M_.nboth   = 1;
M_.nsfwrd   = 5;
M_.nspred   = 3;
M_.ndynamic   = 7;
M_.equations_tags = {
};
M_.exo_names_orig_ord = [1:2];
M_.maximum_lag = 1;
M_.maximum_lead = 1;
M_.maximum_endo_lag = 1;
M_.maximum_endo_lead = 1;
oo_.steady_state = zeros(13, 1);
M_.maximum_exo_lag = 0;
M_.maximum_exo_lead = 0;
oo_.exo_steady_state = zeros(2, 1);
M_.params = NaN(23, 1);
M_.NNZDerivatives = zeros(3, 1);
M_.NNZDerivatives(1) = 49;
M_.NNZDerivatives(2) = 106;
M_.NNZDerivatives(3) = -1;
%
% SHOCKS instructions
%
make_ex_;
M_.exo_det_length = 0;
M_.Sigma_e(1, 1) = 1;
M_.Sigma_e(2, 2) = 1;
M_.sigma_e_is_diagonal = 1;
M_.params( 1 ) = 4.507;
chi = M_.params( 1 );
M_.params( 2 ) = 0.3333333333333333;
eta = M_.params( 2 );
M_.params( 3 ) = 1;
sigma_ = M_.params( 3 );
M_.params( 4 ) = 0.99;
betass = M_.params( 4 );
M_.params( 5 ) = 0.025;
delta = M_.params( 5 );
M_.params( 6 ) = 0.33;
alpha = M_.params( 6 );
M_.params( 7 ) = 58.25;
varphi = M_.params( 7 );
M_.params( 8 ) = 6;
theta = M_.params( 8 );
M_.params( 17 ) = 2.5;
nu = M_.params( 17 );
M_.params( 18 ) = 0.2;
g_over_y = M_.params( 18 );
M_.params( 9 ) = 1.005;
paistar = M_.params( 9 );
M_.params( 10 ) = 1.5;
phi_pai = M_.params( 10 );
M_.params( 11 ) = 0.125;
phi_y = M_.params( 11 );
M_.params( 12 ) = 1;
zbar = M_.params( 12 );
M_.params( 13 ) = 0.8;
rhoz = M_.params( 13 );
M_.params( 14 ) = 0.8;
rhobeta = M_.params( 14 );
M_.params( 15 ) = 0.0025;
sigmabeta = M_.params( 15 );
M_.params( 16 ) = 0.012;
sigmaz = M_.params( 16 );
M_.params( 19 ) = M_.params(9)/M_.params(4);
rbar = M_.params( 19 );
steady;
check(M_,options_,oo_);
options_.order = 2;
var_list_=[];
var_list_ = 'BETA';
var_list_ = char(var_list_, 'C');
var_list_ = char(var_list_, 'I');
var_list_ = char(var_list_, 'K');
var_list_ = char(var_list_, 'N');
var_list_ = char(var_list_, 'PAI');
var_list_ = char(var_list_, 'PSI');
var_list_ = char(var_list_, 'R');
var_list_ = char(var_list_, 'RK');
var_list_ = char(var_list_, 'RTAYLOR');
var_list_ = char(var_list_, 'W');
var_list_ = char(var_list_, 'Y');
info = stoch_simul(var_list_);
save('constant_results.mat', 'oo_', 'M_', 'options_');


disp(['Total computing time : ' dynsec2hms(toc) ]);
diary off
