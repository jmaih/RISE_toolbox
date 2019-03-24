%--- help for ts/dust_up ---
%
%  Sets insignificant digits to zero
% 
%  ::
% 
%    this=dust_up(this)
%    this=dust_up(this,crit)
% 
%  Args:
% 
%     this (rts | ts): time-series object
% 
%     crit (numeric | {sqrt(eps)}): cutoff point
% 
%  Returns:
%     :
% 
%     - **this** [rts|ts]: time-series object
% 
%