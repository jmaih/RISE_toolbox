%  FIRST_ORDER_DERIVATIVES - Compute the first-order partial derivatives of
%  a function or equation. 
% 
%    derivs = first_order_derivatives(ff, wrt, alienList)
% 
%  Input arguments:
%    - ff: A function handle or cell array of function handles representing
%      equations or functions. 
%    - wrt: A cell array specifying variables with respect to which
%      first-order derivatives will be computed. 
%    - alienList (optional): A cell array specifying functions not known to
%      Matlab. 
% 
%  Output:
%    - derivs: A function handle that computes the first-order partial
%      derivatives of the input function(s). 
% 
%  Description:
%  This function calculates the first-order partial derivatives of the
%  provided function(s) with respect to the specified variables in wrt. The
%  derivatives are computed symbolically and returned as a function handle 
%  that can be used to evaluate the derivatives at specific points.
% 
%  Example:
%    % Define equations and variables
%    syms x y z;
%    eq1 = x + y^2;
%    eq2 = z - exp(x);
%    wrt = {'x', 'y', 'z'};
% 
%    % Compute first-order derivatives
%    derivs = first_order_derivatives({eq1, eq2}, wrt);
% 
%    % Evaluate derivatives at a specific point
%    df_dx = derivs(1, 2, 3); % Compute df/dx at x=1, y=2, z=3
% 
%    % The output df_dx contains the value of the derivative.
% 
%  Note: The function computes first-order partial derivatives symbolically
%  and returns a function handle that can be used to evaluate the
%  derivatives at specific points.
%