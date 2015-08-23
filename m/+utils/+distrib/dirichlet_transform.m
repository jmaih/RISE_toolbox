function a=dirichlet_transform(x,a_ii)
% dirichlet_transform -- transforms the parameter of the dirichlet for
% estimation.
%
% Syntax
% -------
% ::
%
%   a=dirichlet_transform(x,a_ii)
%
% Inputs
% -------
%
% - **x** [k-1 x 1 vector]: vector of probabilities excluding the
% "diagonal" element
%
% - **a_ii** [scalar]: un-normalized "weight" of the excluded element.
%
% Outputs
% --------
%
% - **a** [k-1 x 1 vector]: un-normalized "weights" of the off-diagonal
% elements
%
% More About
% ------------
%
% Examples
% ---------
%
% See also: 

% it is assumed that the diagonal element is not part of the vector but can
% be retrieved as
x_i=1-sum(x);

sum_aj=a_ii/x_i-a_ii;

sum_a=a_ii+sum_aj;

a=sum_a*x;

end