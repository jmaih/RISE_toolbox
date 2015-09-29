function aj=dirichlet_transform(x,sum_aij)
% dirichlet_transform -- transforms the parameter of the dirichlet for
% estimation.
%
% Syntax
% -------
% ::
%
%   aj=dirichlet_transform(x,sum_aij)
%
% Inputs
% -------
%
% - **x** [k-1 x 1 vector]: vector of probabilities excluding the
% "diagonal" element
%
% - **sum_aij** [scalar]: sum of the weights including the "diagonal" element.
%
% Outputs
% --------
%
% - **aj** [k-1 x 1 vector]: un-normalized "weights" of the off-diagonal
% elements
%
% More About
% ------------
%
% Examples
% ---------
%
% See also: 

aj=sum_aij*x;

end