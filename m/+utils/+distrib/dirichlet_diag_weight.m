function a_ii=dirichlet_diag_weight(h,p_ii_min,a_ij_max)
% dirichlet_diag_weight -- computes the un-normalized "weight" of the
% missing diagonal element in the transformed/inverse dirichlet
%
% Syntax
% -------
% ::
%
%   a_ii=dirichlet_diag_weight(h)
%
%   a_ii=dirichlet_diag_weight(h,p_ii)
%
%   a_ii=dirichlet_diag_weight(h,p_ii,a_ij_max)
%
% Inputs
% -------
%
% - **h** [integer]: number of regimes
%
% - **p_ii_min** [scalar|vector|{0.05}]: minimum probability of staying in
% regime "i".
%
% - **a_ij_max** [scalar|vector|{1}]: maximum value for the un-normalized
% "weight" of the off-diagonal terms 
%
% Outputs
% --------
%
% - **a_ii** [numeric]: un-normalized "weight" of the reference entry
% (diagonal term) 
%
% More About
% ------------
%
% - Under estimation then, the transformed parameters will be
% xj=aj/(a_ii+sum(aj)), where aj is in [0,a_ij_max]
%
% Examples
% ---------
%
% See also: 

if nargin<3
    a_ij_max=[];
    if nargin<2
        p_ii_min=[];
    end
end
if isempty(a_ij_max)
    a_ij_max=1;
end
if isempty(p_ii_min)
    p_ii_min=0.05;
end
a_ii=p_ii_min./(1-p_ii_min).*(h-1).*a_ij_max;
end