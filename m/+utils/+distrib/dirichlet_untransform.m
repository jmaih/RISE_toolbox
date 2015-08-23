function x=dirichlet_untransform(a,a_ii)
% dirichlet_untransform -- sets the transformed dirichlet parameters back
% to probabilities
%
% Syntax
% -------
% ::
%
%   x=dirichlet_untransform(a,a_ii)
%
% Inputs
% -------
%
% - **a** [k-1 x 1 vector]: un-normalized "weights" of the off-diagonal
% elements
%
% - **a_ii** [scalar]: un-normalized "weight" of the excluded element.
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

x=a/(a_ii+sum(a));
end