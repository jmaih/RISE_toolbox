%--- help for ts/npdecomp ---
%
%  Non-parametric decomposition into trend, seasonal and irregular components
% 
%  ::
% 
%    out=npdecomp(y)
%    out=npdecomp(y,doLog)
% 
%  Args:
% 
%     y (ts): time series to decompose
% 
%     doLog (true | {false}): if log, do a multiplicative decomposition
%       otherwise the decomposition is additive
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
%     - :func:`pdecomp <ts.pdecomp>`
% 
%