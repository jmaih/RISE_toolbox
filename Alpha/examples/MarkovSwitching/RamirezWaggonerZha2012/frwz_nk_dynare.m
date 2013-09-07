%
% Status : main Dynare file 
%
% Warning : this file is generated automatically by Dynare
%           from model file (.mod)

tic;
global M_ oo_ options_ ys0_ ex0_ estimation_info
options_ = [];
M_.fname = 'frwz_nk_dynare';
%
% Some global variables initialization
%
global_initialization;
diary off;
logname_ = 'frwz_nk_dynare.log';
if exist(logname_, 'file')
    delete(logname_)
end
diary(logname_)
M_.exo_names = 'EPS_R';
M_.exo_names_tex = 'EPS\_R';
M_.endo_names = 'PAI';
M_.endo_names_tex = 'PAI';
M_.endo_names = char(M_.endo_names, 'Y');
M_.endo_names_tex = char(M_.endo_names_tex, 'Y');
M_.endo_names = char(M_.endo_names, 'R');
M_.endo_names_tex = char(M_.endo_names_tex, 'R');
M_.param_names = 'betta';
M_.param_names_tex = 'betta';
M_.param_names = char(M_.param_names, 'eta');
M_.param_names_tex = char(M_.param_names_tex, 'eta');
M_.param_names = char(M_.param_names, 'kappa');
M_.param_names_tex = char(M_.param_names_tex, 'kappa');
M_.param_names = char(M_.param_names, 'mu');
M_.param_names_tex = char(M_.param_names_tex, 'mu');
M_.param_names = char(M_.param_names, 'mu_bar');
M_.param_names_tex = char(M_.param_names_tex, 'mu\_bar');
M_.param_names = char(M_.param_names, 'psi');
M_.param_names_tex = char(M_.param_names_tex, 'psi');
M_.param_names = char(M_.param_names, 'rhor');
M_.param_names_tex = char(M_.param_names_tex, 'rhor');
M_.param_names = char(M_.param_names, 'sigr');
M_.param_names_tex = char(M_.param_names_tex, 'sigr');
M_.exo_det_nbr = 0;
M_.exo_nbr = 1;
M_.endo_nbr = 3;
M_.param_nbr = 8;
M_.orig_endo_nbr = 3;
M_.aux_vars = [];
M_.Sigma_e = zeros(1, 1);
M_.H = 0;
options_.block=0;
options_.bytecode=0;
options_.use_dll=0;
erase_compiled_function('frwz_nk_dynare_static');
erase_compiled_function('frwz_nk_dynare_dynamic');
M_.lead_lag_incidence = [
 0 2 5;
 0 3 6;
 1 4 0;]';
M_.nstatic = 0;
M_.nfwrd   = 2;
M_.npred   = 1;
M_.nboth   = 0;
M_.nsfwrd   = 2;
M_.nspred   = 1;
M_.ndynamic   = 3;
M_.equations_tags = {
};
M_.exo_names_orig_ord = [1:1];
M_.maximum_lag = 1;
M_.maximum_lead = 1;
M_.maximum_endo_lag = 1;
M_.maximum_endo_lead = 1;
oo_.steady_state = zeros(3, 1);
M_.maximum_exo_lag = 0;
M_.maximum_exo_lead = 0;
oo_.exo_steady_state = zeros(1, 1);
M_.params = NaN(8, 1);
M_.NNZDerivatives = zeros(3, 1);
M_.NNZDerivatives(1) = 12;
M_.NNZDerivatives(2) = 38;
M_.NNZDerivatives(3) = -1;
a_tp_1_2=1-.9; 
a_tp_2_1=1-.9;
M_.params( 1 ) = .99;
betta = M_.params( 1 );
M_.params( 3 ) = 161;
kappa = M_.params( 3 );
M_.params( 2 ) = 10;
eta = M_.params( 2 );
M_.params( 7 ) = .8;
rhor = M_.params( 7 );
M_.params( 8 ) = 0.0025;
sigr = M_.params( 8 );
M_.params( 5 ) = 0.02;
mu_bar = M_.params( 5 );
M_.params( 4 ) = M_.params(5);
mu = M_.params( 4 );
M_.params( 6 ) = 3.1;
psi = M_.params( 6 );
steady;
check(M_,options_,oo_);
%
% SHOCKS instructions
%
make_ex_;
M_.exo_det_length = 0;
M_.Sigma_e(1, 1) = 1;
M_.sigma_e_is_diagonal = 1;
options_.irf = 40;
options_.order = 2;
var_list_=[];
info = stoch_simul(var_list_);
save('frwz_nk_dynare_results.mat', 'oo_', 'M_', 'options_');


disp(['Total computing time : ' dynsec2hms(toc) ]);
diary off
