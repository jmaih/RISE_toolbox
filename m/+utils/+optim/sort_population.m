function [pop,order]=sort_population(pop,iteration)
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

if nargin<2
    iteration=[];
end
if isempty(iteration)
    % Deb sorting
    violstrength=[pop.violstrength];
    fs=[pop.f];
    order=[];
    uv=unique(violstrength);
    for iv=1:numel(uv)
        locs=find(violstrength==uv(iv));
        [~,tags]=sort(fs(locs));
        order=[order,locs(tags)]; %#ok<AGROW>
    end
else
    alpha_=1000;
    % apply dynamic penalty and then sort
    nvec=numel(pop);
    tt=nan(1,nvec);
    for ii=1:nvec
        tt(ii)=pop(ii).f+utils.optim.dynamic_penalty(pop(ii).viol,iteration,alpha_);
    end
    [~,order]=sort(tt);
end
pop=pop(order);
end
