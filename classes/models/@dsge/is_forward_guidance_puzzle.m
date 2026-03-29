%--- help for dsge/is_forward_guidance_puzzle ---
%
%  is_forward_guidance_puzzle: Checks whether a DSGE model exhibits a
%  forward-guidance puzzle. 
% 
%  The function can handle multiple models simultaneously, and each model
%  can have multiple parameterizations. Under regime switching, it is not
%  possible to characterize the forward guidance as in constant-parameter
%  models because the impact of shocks differs across regimes.
% 
%  Syntax::
% 
%    o=is_forward_guidance_puzzle(m)
% 
%    o=is_forward_guidance_puzzle(m,errIfMultReg)
% 
%  Inputs:
% 
%    - m: DSGE model or an array of DSGE models to analyze.
% 
%    - errIfMultReg [{true}|false]: If ``True``, raise an error when handling
%      regime switching models. If ``False``, assume shock impacts are the same
%      across regimes.
% 
%  Output:
% 
%    - o: A logical array indicating whether each model exhibits a
%      forward-guidance puzzle for each parameterization and regime.
% 
%    - maxEig: maximum eigenvalue criterion
% 
%  The function checks whether a given DSGE model or a list of DSGE models
%  exhibit a forward-guidance puzzle. The forward-guidance puzzle occurs
%  when the model's dynamics lead to explosive responses to future shocks.
% 
%  This function returns a logical array where each element corresponds to a
%  specific model parameterization and regime. A value of ``True`` indicates
%  that a forward-guidance puzzle is present, and ``False`` indicates no
%  forward-guidance puzzle.
% 
%  .. note::
% 
%    - When analyzing a single model, the result is a logical array.
%    - When analyzing multiple models, the result is a cell array of logical
%      arrays, with one element for each model. 
% 
%  .. warning::
% 
%    When dealing with regime-switching models (models with multiple regimes),
%    setting ``errIfMultReg=True`` will raise an error because it's not
%    possible to characterize forward guidance in such models due to differing
%    shock impacts across regimes.
%