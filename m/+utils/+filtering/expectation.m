function E=expectation(probs,vals,straight)
% EXPECTATION -- computes the expectation of variables from multiple regimes
%
% ::
%
%
%   E=expectation(probs,vals)
%
% Args:
%
%    - **probs** [k x T matrix]: matrix of probabilities
%
%    - **vals** [1 x k cell]: matrix of values for which we want to take the
%    expectation
%
%    - **straight** [true|{false}]: decides whether permutations should be
%    taken or not.
%
% Returns:
%    :
%
%    - **E** [1 x T vector]: expecation
%
% Note:
%
% Example:
%
%    See also:

if nargin<3
    % legacy for DSGE... could change in the future
    straight=false;
end

E=0;
for istate=1:numel(vals)
    if straight
        % VARs
        E=E+bsxfun(@times,probs(istate,:),vals{istate});
    else
        % DSGEs
        % E=E+squeeze(bsxfun(@times,permute(probs(istate,:),[3,1,2]),vals{istate}));
        E=E+bsxfun(@times,permute(probs(istate,:),[3,1,2]),vals{istate});
    end
end
% E=squeeze(sum(bsxfun(@times,permute(probs,[3,1,2]),vals),2));
end

