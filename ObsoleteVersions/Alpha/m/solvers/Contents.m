%%%%%%%%%%%%%%%%%%%%   path: m\solvers   %%%%%%%%%%%%%%%%%%%%
%
%   <a href="matlab:help m\solvers\collapse_array">m\solvers\collapse_array</a>                               - takes as input a multidimensional array and outputs a cell array
%   <a href="matlab:help m\solvers\computational_savings">m\solvers\computational_savings</a>                        - this function separates static variables from dynamic ones. It places all
%   <a href="matlab:help m\solvers\dsge_lc_solve">m\solvers\dsge_lc_solve</a>                                - Reference: Debortoli, Maih, Nunes (2010): "Loose Commitment in Medium
%   <a href="matlab:help m\solvers\dsge_solve_aim">m\solvers\dsge_solve_aim</a>                               - ags=1; % no of lags and leads
%   <a href="matlab:help m\solvers\dsge_solve_gensys">m\solvers\dsge_solve_gensys</a>                            - this function solves the rational expectations model
%   <a href="matlab:help m\solvers\dsge_solve_klein">m\solvers\dsge_solve_klein</a>                             - this function solves the rational expectations model
%   <a href="matlab:help m\solvers\expand_array">m\solvers\expand_array</a>                                 - case 1	% A0, Aminus or Aplus
%   <a href="matlab:help m\solvers\fix_point_iterator">m\solvers\fix_point_iterator</a>                           - this function solves for a fix point. Inputs are:
%   <a href="matlab:help m\solvers\gensys">m\solvers\gensys</a>                                       - function [G1,C,impact,fmat,fwt,ywt,gev,eu,loose]=gensys(g0,g1,c,psi,pi,div)
%   <a href="matlab:help m\solvers\get_default_optimization_option">m\solvers\get_default_optimization_option</a>              - lc_reconvexify:
%   m\solvers\is_eigenvalue_solver_candidate               - (No help available)
%   <a href="matlab:help m\solvers\loose_commitment_solver">m\solvers\loose_commitment_solver</a>                      - Reference: Debortoli, Maih, Nunes (2010): "Loose Commitment in Medium
%   <a href="matlab:help m\solvers\loose_commitment_solver_fix_point_unfinished">m\solvers\loose_commitment_solver_fix_point_unfinished</a> - Reference: Debortoli, Maih, Nunes (2010): "Loose Commitment in Medium
%   <a href="matlab:help m\solvers\loose_commitment_to_markov_switching">m\solvers\loose_commitment_to_markov_switching</a>         - this function puts the loose commitment solution into a markov switching
%   m\solvers\markov_switching_dsge_objective              - (No help available)
%   <a href="matlab:help m\solvers\markov_switching_dsge_stack">m\solvers\markov_switching_dsge_stack</a>                  - C_st,... % endo_nbr x h matrix of constant
%   <a href="matlab:help m\solvers\msre_aim">m\solvers\msre_aim</a>                                     - solve for T only
%   <a href="matlab:help m\solvers\msre_gensys">m\solvers\msre_gensys</a>                                  - solve for T only
%   <a href="matlab:help m\solvers\msre_klein">m\solvers\msre_klein</a>                                   - solve for T only
%   <a href="matlab:help m\solvers\msre_matrix_times_vector">m\solvers\msre_matrix_times_vector</a>                     - the old version is faster and solves but solves the problem
%   <a href="matlab:help m\solvers\msre_solve">m\solvers\msre_solve</a>                                   - This procedure assumes the steady state has been solved and that apart
%   <a href="matlab:help m\solvers\qzdiv">m\solvers\qzdiv</a>                                        - function [A,B,Q,Z] = qzdiv(stake,A,B,Q,Z)
%   <a href="matlab:help m\solvers\qzswitch">m\solvers\qzswitch</a>                                     - function [A,B,Q,Z] = qzswitch(i,A,B,Q,Z)
%   <a href="matlab:help m\solvers\schur_solver">m\solvers\schur_solver</a>                                 - T,S,Q,Z] = qz(F,G,'real');%complex
%   <a href="matlab:help m\solvers\solve_steady_state">m\solvers\solve_steady_state</a>                           - %         ys0=0*ys0;
%   <a href="matlab:help m\solvers\transpose_free_quasi_minimum_residual">m\solvers\transpose_free_quasi_minimum_residual</a>        - A,... % coefficient matrix
%   <a href="matlab:help m\solvers\unearth_frwrd_matrix">m\solvers\unearth_frwrd_matrix</a>                         - this function extracts Aplus from Gplus and is used for solving models
