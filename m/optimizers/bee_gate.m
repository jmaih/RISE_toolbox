%  `BEE_GATE` attempts to find the global minimum of a constrained function using the BEE algorithm.
% 
%    [x, f, exitflag, H, obj] = bee_gate(Objective, x0, lb, ub, options, varargin)
% 
%    Inputs:
%    - `Objective`: Function handle to the objective function.
%    - `x0`: Initial guess for the optimization.
%    - `lb`: Lower bounds for the variables.
%    - `ub`: Upper bounds for the variables.
%    - `options`: Options for the BEE algorithm. (Optional)
%    - `varargin`: Additional arguments for the objective function.
% 
%    Outputs:
%    - `x`: Optimized solution.
%    - `f`: Objective function value at the optimized solution.
%    - `exitflag`: Exit flag, indicating the reason for termination.
%    - `H`: Estimate of the Hessian at the optimized solution.
%    - `obj`: Optimized BEE object with information about the optimization process.
% 
%    Options Structure Fields:
%    - `MaxNodes`: Number of different elements in the group sharing information. (Default: 20)
%    - `MaxIter`: Maximum number of iterations. (Default: 1000)
%    - `MaxTime`: Time budget in seconds. (Default: 3600)
%    - `MaxFunEvals`: Maximum number of function evaluations. (Default: inf)
%    - `rand_seed`: Seed number for random draws.
% 
%    Notes:
%    - Optimization stops when one of the following happens:
%      1. The number of iterations exceeds MaxIter.
%      2. The number of function counts exceeds MaxFunEvals.
%      3. The time elapsed exceeds MaxTime.
%      4. The user writes anything in and saves the automatically generated file called "ManualStopping.txt".
% 
%    Examples:
%      FUN can be specified using `@`:
%        x = bee_gate(@myfunc, ...)
%      In this case, `F = myfunc(X)` returns the scalar function value `F` of the MYFUNC function evaluated at `X`.
% 
%      FUN can also be an anonymous function:
%        x = bee_gate(@(x) 3*sin(x(1))+exp(x(2)), [1;1], [], [], [], [], [0 0])
%      returns `X = [0;0]`.
% 
%      FUN = inline('sum(x.^2)'); n = 100;
%      lb = -20*ones(n,1); ub = -lb; x0 = lb + (ub - lb) .* rand(n,1);
%      options = struct('MaxNodes', 20, 'MaxIter', 1000, 'MaxTime', 60, 'MaxFunEvals', inf);
%      x = bee_gate(@(x) FUN(x), x0, lb, ub, options);
% 
%    See also: BEE
% 
%    Copyright 2011 Junior Maih (junior.maih@gmail.com).
%    $Revision: 7 $  $Date: 2011/05/26 11:23 $
%