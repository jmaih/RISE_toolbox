function covariance_tune=retune_coefficient(covariance_tune,target_range,acceptance_rate)
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

if acceptance_rate<target_range(1)||acceptance_rate>target_range(2)
    % increase c if the acceptance rate is high
    % decrease it otherwise
    target_rate=mean(target_range);
    covariance_tune=covariance_tune*acceptance_rate/target_rate;
end
end