%--- help for ts/chebyshev_box ---
%
%  Constructs chebyshev boxes for multivariate-multiperiods densities
% 
%  ::
% 
%    mvcb=chebyshev_box(this,gam)
% 
%  Args:
% 
%     this (ts | rts) : time series with many pages (number of simulations)
%       and potentially many columns (number of variables)
% 
%     gam (scalar | vector) : percentile(s)
% 
%  Returns:
%     :
% 
%     - **mvcb** [struct] : structure with time series in its fields. Each
%       field represents a particular variable
% 
%     - **gam_** [scalar|vector] : sorted percentile(s)
% 
%     - **my** [ts] : mean across simulations
% 
%