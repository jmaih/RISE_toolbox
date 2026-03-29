%--- help for dsge/historical_decomposition_switch ---
%
%  historical_decomposition Computes historical decomposition of a
%  regime-switching DSGE model 
% 
%  Syntax::
% 
%    mycontrib=historical_decomposition_switch(m)
% 
%    mycontrib=historical_decomposition_switch(m,groups)
% 
%    mycontrib=historical_decomposition_switch(m,groups,varargin)
% 
%  Inputs:
% 
%  - m : [rise|dsge] model(s) for which to compute the historical 
%    decomposition. m could be a vector of models. The computation of the
%    historical decomposition requires a computation of smoothed history.
%    Hence either m should contain all the information needed for that
%    purpose or the information should be passed along through varargin below
% 
%  - groups : [structure|cell array |{empty}] grouping of shocks in the decomposition.
%    By default, the shocks are not grouped. The syntax is of the form
%    {group1,{v11,v12,...},...,groupn,{vn1,vn2,...}}. The shocks that are
%    not listed are put in a special group called "others". The "others"
%    group does not include the effect of initial conditions.
%    e.g. p=struct(); p.demand={'Ey','Er'}; p.supply={'Ep'};
%    e.g. p={'demand',{'Ey','Er'},'supply',{'Ep'}};
% 
%  - varargin : additional information needed for the computation of the
%    smoothed quantities through filtration.
% 
%  Outputs:
% 
%  - mycontrib : [struct|cell array] structure or cell array of structures
%    with the contributions in each model and each solution. The
%    decompositions are given in terms of:
% 
%    - the exogenous variables
%    - **init** : the effect of initial conditions, which includes the
%      steady state!!!
% 
%  Remarks:
% 
%  - **N.B** : For a model with switching, RISE adds the contribution of the
%    switching process
% 
%  - **N.B** : if m is a vector of models, then each model should return a
%    unique solution (implying a unique filtration), else the concatenation
%    of decompositions will fail. In that case it is better to run one model
%    at a time.
% 
%  - **N.B** : When a model has multiple solutions, there is no guarantee
%    that every solution will have succeed at filtration. For this reason,
%    the result of each decomposition is stored in a separate cell. Empty
%    cells represent the solutions for which the filtration could not be
%    obtained.
% 
%  - **N.B** : For the variables declared as "log_variables", the
%     decomposition is given in terms of transformed variables and not in
%     terms of their level or original units. 
% 
%  Examples:
% 
%  See also: dsge/historical_decomposition
%