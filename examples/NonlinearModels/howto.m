%% housekeeping
close all
clear
clc
%% add the necessary paths
rise_startup()

%% we read the model
m=dsge('fs2000');
%% we solve, passing a function that solves the steady state analytically

mwssf_a=solve(m,'steady_state_file','fs2000_steadystate');

%% we solve passing a function that gives initial values of the steady state. 
% this is the counterpart to dynare's initval

mwssf_b=solve(m,'steady_state_file','fs2000_steadystate_initval');

%% compare the solutions
mwssf_a.print_solution
mwssf_b.print_solution

%% now we repeat the same exercise but with the steady state model inside the model file
% since we would like to do the two cases above in an elegant way, we
% introduce a switch called 'approximate'. You can call it what you want
% but it has to be issued while 'rising' the model. Not after.

% to continue in our elegance we initialize a vector of models
mv=rise.empty(0);
approximate=[0,1];
for approx=1:numel(approximate)
    approx_flag=struct('approximate',approximate(approx));
    % when approximate==0, we read the exact steady state solution
    % when approximate==1, we read an approximation of the steady state
    % solution. in both cases, the computed values will be used as an
    % initial guess for the computation of the steady state, unless we add
    % the attribute 'imposed' to 'steady_state_model', in which case rise
    % does not check that the initial guess actually solves the steady
    % state.
    mv(approx,1)=rise('fs2000_b','rise_flags',approx_flag);
end

%% we solve both models simultaneously
mv=solve(mv);

%% we print the solution
print_solution(mv)

%% a word of warning
% you might be tempted to or make the mistake of having both a
% steady_state_model block and a steady state file. RISE will use the
% steady state file and ignore the steady_state_model block!!! If you don't
% like this behavior, shoot an email to junior.maih AT gmail.com and let me
% know about your arguments for doing it differently.