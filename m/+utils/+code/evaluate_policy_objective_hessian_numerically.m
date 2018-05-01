function H=evaluate_policy_objective_hessian_numerically(funcs,y,x,ss,param,sparam,def,s0,s1)
% evaluate_policy_objective_hessian_numerically - numerical evaluation of the hessian of the policy objective
%
% ::
%
%
%   H=evaluate_policy_objective_hessian_numerically(funcs,y,x,ss,param,sparam,def,s0,s1)
%
% Args:
%
%    - **funcs** [fhandle|cell array]: function or functions to be
%      differentiated
%
%    - **y** [vector]: values of endogenous variables
%
%    - **x** [vector]: values of exogenous variables
%
%    - **ss** [vector]: steady state
%
%    - **param** [vector]: parameter vector
%
%    - **sparam** [vector]: vector of parameters appearing with a lead
%
%    - **def** [vector]: values of definitions
%
%    - **s0** [scalar]: state today
%
%    - **s1** [scalar]: state tomorrow
%
% Returns:
%    :
%
%    - **H** [matrix]: Numerical Hessian of **funcs** at [y]
%
% Note:
%
%    It is assumed that the inputs are y,x,ss,param,sparam,def,s0,s1 as
%    ordered in parser.input_list()
%
% Example:
%
%    See also:

if iscell(funcs)
    funcs=funcs{1};
end

if ~isa(funcs,'function_handle')
    error('objective should be a function handle')
end

H=utils.numdiff.hessian(@newobjective,y);

    function f=newobjective(yy)
        f=funcs(yy,x,ss,param,sparam,def,s0,s1);
    end
end