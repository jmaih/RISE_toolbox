function [a,P,eta]=squasher(a,P,m)
% squasher - separates filters into shocks and filters
%
% Syntax
% -------
% ::
%
%   [a,P,eta]=squasher(a,P,m)
%
% Inputs
% -------
%
% - **a** [3-D array]: filtered, updated or smoothed variables INcluding
%   shocks
%
% - **P** [3-D array]: covariance matrix for filtered, updated or smoothed
%   variables INcluding shocks
%
% - **m** [scalar]: number of endogenous variables
%
% Outputs
% --------
%
% - **a** [3-D array]: filtered, updated or smoothed variables EXcluding
%   shocks
%
% - **P** [3-D array]: covariance matrix for filtered, updated or smoothed
%   variables EXcluding shocks
%
% - **eta** [3-D array]: filtered, updated or smoothed shocks
%
% More About
% ------------
%
% Examples
% ---------
%
% See also: 

eta=a(m+1:end,1,:);
a=a(1:m,:,:); % second dimension is the real-time forecasting steps
P=P(1:m,1:m,:);

end