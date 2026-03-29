%  svar : constructor for structural VAR models
% 
%  Syntax :
% 
%  ::
% 
%        mdl=svar(varlist)
% 
%        mdl=svar(varlist,exog)
% 
%        mdl=svar(varlist,exog,nlags)
% 
%        mdl=svar(varlist,exog,nlags,constant)
% 
%        mdl=svar(varlist,exog,nlags,constant,markov_chains)
% 
%  Args :  
% 
%  - **varlist** (cellstr): List of endogenous variables 
% 
%  - **exog** (empty|cellstr): List of exogenous variables
%    excluding the constant term    
% 
%  - **nlags** (empty|integer|{4}): Number of lags in the VAR 
% 
%  - **constant** (empty|{true}|false): boolean for including a
%    constant term in the VAR or not    
% 
%  - **markov_chains** (empty|struct): structure containing the following fields:
% 
%        - **'name'**: name of the markov process
%        - **'number_of_states'**: number of states
%        - **'controlled_parameters'**: list of controlled
%          parameters
%        - **'endogenous_probabilities'**: definitions of the
%          endogenous probabilities
%        - **'probability_parameters'**: parameters entering the
%          time-varying transition probabilities 
% 
%  Returns :  
% 
%  - **mdl** (svar): Constructed SVAR model 
% 
%  .. note ::
% 
%    For the controlled parameters
% 
%    - The dynamic parameters of a svar model are denoted by "a"
%    - The static parameters or coefficients on exogenous by "c"
%    - for some parameter bd(row,col), e.g. a2(2,3)
% 
%        - d (empty|integer) : when empty, all lags are swept
%          through
%        - row (integer|string) : denotes the row, which can be
%          an integer of a string corresponding to the name of
%          an endogenous variable
%        - col (integer|string) : denotes the row, which can be
%          an integer of a string corresponding to the name of
%          an endogenous variable
%
%    Documentation for svar
%       doc svar
%
%