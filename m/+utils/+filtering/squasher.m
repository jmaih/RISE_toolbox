%  INTERNAL FUNCTION: separates filters into shocks and filters
% 
%  ::
% 
%    [a,P,eta]=squasher(a,P,m)
% 
%  Args:
% 
%     - **a** [3-D array]: filtered, updated or smoothed variables INcluding
%       shocks
%     - **P** [3-D array]: covariance matrix for filtered, updated or smoothed
%       variables INcluding shocks
%     - **m** [scalar]: number of endogenous variables
% 
%  Returns:
%     :
% 
%     - **a** [3-D array]: filtered, updated or smoothed variables EXcluding
%       shocks
%     - **P** [3-D array]: covariance matrix for filtered, updated or smoothed
%       variables EXcluding shocks
%     - **eta** [3-D array]: filtered, updated or smoothed shocks
% 
%