%--- help for ts/pdecomp ---
%
%  Parametric decomposition into trend, seasonal and irregular components
% 
%  ::
% 
%    out=pdecomp(y)
%    out=pdecomp(y,doLog)
%    out=pdecomp(y,doLog,dorder)
% 
%  Args:
% 
%     y (ts): time series to decompose
%     doLog (true | {false}): if log, do a multiplicative decomposition
%       otherwise the decomposition is additive
%     dorder (integer | {2}): detrending order
% 
%  Returns:
%     :
% 
%     - **out** [struct] :
% 
%       - **trend** [ts] : estimated trend
%       - **sc**    [ts] : estimated seasonal component
%       - **sa**    [ts] : seasonally adjusted data
%       - **ic**    [ts] : estimated irregular component
% 
%  Note:
% 
%     If there are many variables and the variables are named, the first level
%     of the structure will be the names of the different variables.
% 
%  See also:
%     - :func:`npdecomp <ts.npdecomp>`
% 
%