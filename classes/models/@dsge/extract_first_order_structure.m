%--- help for dsge/extract_first_order_structure ---
%
%  EXTRACT_FIRST_ORDER_STRUCTURE Extracts the first-order structure of a DSGE model.
% 
%  Syntax:
% 
%    [Aplus, A0, Aminus, B, Q, ss, growth] = extract_first_order_structure(obj)
%    [Aplus, A0, Aminus, B, Q, ss, growth] = extract_first_order_structure(obj, sm)
% 
%  Inputs:
% 
%    - obj [rise|dsge]: Model object.
%    - sm [struct]: Derivatives computed by RISE (optional).
% 
%  Outputs:
% 
%    - Aplus [h x h cell array]: Coefficients on the forward-looking terms
%      multiplied by the transition matrix. 
%    - A0 [1 x h cell array]: Coefficients on the contemporaneous terms.
%    - Aminus [1 x h cell array]: Coefficients on the backward-looking terms.
%    - B [1 x h cell array]: Coefficients on shocks.
%    - Q [h x h matrix]: Transition matrix (qij) with i=today and j=tomorrow.
%    - ss [n x h matrix]: Steady states.
%    - growth [n x h matrix]: Balanced growth.
% 
%  Notes:
% 
%    - The steady state is returned instead of the constant term. This
%      assumes that the model is at least conditionally stationary. The
%      constant term can be recovered as b = -(Aplus{i,i}/Q(i,i) + A0{i} +
%      Aminus{i}) * ss(:, i).    
%    - In case of a model with multiple parameterizations, all the outputs
%      are cell arrays, with each cell representing one parameterization. 
% 
%  Description:
% 
%    The function extracts the first-order structure of a DSGE model,
%    allowing the user to use their own algorithms externally. It returns
%    the coefficients on various terms in the linearized model equations, as
%    well as the transition matrix, steady states, and balanced growth rates.
% 
%  Example:
%    % Extract first-order structure from a DSGE model
%    model = rise('my_model_file');
%    [Aplus, A0, Aminus, B, Q, ss, growth] = extract_first_order_structure(model);
% 
%  See also:
%    solve
%