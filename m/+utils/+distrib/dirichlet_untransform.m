function x=dirichlet_untransform(aj,sum_aij)
% dirichlet_untransform -- sets the transformed dirichlet parameters back
% to probabilities
%
% Syntax
% -------
% ::
%
%   x=dirichlet_untransform(aj,sum_aij)
%
% Inputs
% -------
%
% - **aj** [k-1 x 1 vector]: un-normalized "weights" of the off-diagonal
% elements
%
% - **sum_aij** [scalar]: sum of "weights" including the "diagonal" element.
%
% Outputs
% --------
%
% - **x** [k-1 x 1 vector]: vector of probabilities excluding the
% "diagonal" element
%
% More About
% ------------
%
% Examples
% ---------
%
% See also: 

x=aj/sum_aij;

end