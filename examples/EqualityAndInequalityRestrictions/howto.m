%% inequality restrictions
% use the parameter restrictions block
% a >= b
% cos(tan(b))<=log(c)

%% equality restrictions
% One can use the parameter restrictions block as the interface
% to equality restrictions. Those restrictions can be nonlinear
% a = 3*b+4*c
% a = tan(cos(log(sin(b+c))))
% in the example above, if a is not estimated while b and c can be
%
% An alternative and more general solution is to use a steady state file or
% the steady_state_model block.

%% inequality restrictions as equality restrictions
% suppose we have the following inequality restriction
% a>=b
% this implies that there is a c>=0 such that a=b+c. One can exploit this
% insight to transform inequality restrictions into equality restrictions.
% We would need to define b and c as estimated while a will not appear in
% the parameterization block. The implied equality restriction will then be
% written either in the parameter_restrictions block (when this is formally
% implemented) or in a steady state file

%% fmincon and estimation with inequality restrictions
% lately, my experience with fmincon is that it does not always handle
% inequality restrictions well. Estimation may break down completely while
% fmincon tries to compute the derivatives of the inequality restrictions
% at a particular point. Unfortunately the function doing that step is
% p-coded and it is difficult to implement a workaround within fmincon
% itself (The Mathworks should be able to help on that...)
% in the models we estimate, this is sometimes the case for transition
% probabilities. Fmincon is known to work well when problems are smooth and
% differentiable. So in order to avoid the breakdown of fmincon, one has to
% be smart and RISE is:
% 1- use priors to smooth the shape of the objective function
% 2- use block_wise optimization
% 3- use a different optimization algorithm like bee_gate and its brothers,
% which I still need to push into the current version of RISE
%{
myblocks={
    {'beta_trans','kappa','rho_d','rho_s','rho_o','Lambda','Gamma',...
    'sig_d(vol,1)','sig_d(vol,2)','sig_s(vol,1)','sig_s(vol,2)',...
    'sig_r(vol,1)','sig_r(vol,2)','sig_o(oil,1)','sig_o(oil,2)','phi_pi',...
    'phi_y','rho_r'}
    {'vol_tp_1_2','vol_tp_2_1','oil_tp_1_2','oil_tp_2_1'}
};
m = estimate(m, 'optimizer', 'fmincon',...
    'estim_blocks',myblocks);
%}