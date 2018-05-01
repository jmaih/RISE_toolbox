function [x,theta_low,theta_high]=theta_to_x(theta,theta_low,theta_high,xlow,xhigh)
% H1 line
%
% ::
%
%
% Args:
%
% Returns:
%    :
%
% Note:
%
% Example:
%
%    See also:

if nargin<5
    xhigh=1;
    if nargin<4
        xlow=0;
        if nargin<3
            theta_high=[];
            if nargin<2
                theta_low=[];
                if nargin<1
                    error([mfilename,':: wrong number of arguments'])
                end
            end
        end
    end
end

% normalize
if isempty(theta_low)
    theta_low=min(theta,[],2);
end
if isempty(theta_high)
    theta_high=max(theta,[],2);
end
x=xlow+(bsxfun(@rdivide,bsxfun(@minus,theta,theta_low),theta_high-theta_low))*(xhigh-xlow);
