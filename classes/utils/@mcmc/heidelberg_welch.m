function [p]=heidelberg_welch(obj,alpha)

if nargin<2
    
    alpha=[];
    
end

if isempty(alpha)
    
    alpha=0.05;
    
end

error('Implementation unfinished')

% First part
%-----------
% 2. Calculate the test statistic on the whole chain. Accept or reject null
% hypothesis that the chain is from a stationary distribution
%
% 3. If null is rejected, discard the first 10% of the chain. Calculate the
% test statistic and accept or reject null... discard 20%, 30%, 40% 50% if
% necessary to find the stationary distribution
%

% Second part
%------------
% 2. Calculate the test statistic on the whole chain. Accept or reject null
% hypothesis that the chain is from a stationary distribution
%
% 3. If null is rejected, discard the first 10% of the chain. Calculate the
% test statistic and accept or reject null... discard 20%, 30%, 40% 50% if
% necessary to find the stationary distribution
%

end