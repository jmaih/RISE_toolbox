function Island=recenter(Island,lb,ub,flag)
% RECENTER - reposition draws without their boundaries
%
% ::
%
%
%   Island=RECENTER(Island,lb,ub)
%
%   Island=RECENTER(Island,lb,ub,flag)
%
% Args:
%
%    - **Island** [vector]: draw to recenter
%
%    - **lb** [vector]: lower bound of search space
%
%    - **ub** [vector]: upper bound of search space
%
%    - **flag** [{1}|2|3]:
%      - **1**: set offending parameters to the bounds they exceed
%      - **2**: set offending parameters to min(ub,2*lb-param) if they are
%      below the lower bound and to max(lb,2*ub-param) if they exceed the
%      upper bound.
%      - **3**: set offending parameters to random draws withing the bounds
%
% Returns:
%    :
%
%    - **Island** [vector]: recentered parameter vector
%
% Note:
%
% Example:
%
%    See also:

if nargin<4
    flag=1;
end
pop=size(Island,2);
lb=lb(:,ones(pop,1));
ub=ub(:,ones(pop,1));
negs=Island<lb;
pos=Island>ub;
switch flag
    case 1
        Island(negs)=lb(negs);
        Island(pos)=ub(pos);
    case 2
        Island(negs)=min(ub(negs),2*lb(negs)-Island(negs));
        Island(pos)=max(lb(pos),2*ub(pos)-Island(pos));
    case 3
        if any(negs)
            Island(negs)=lb(negs)+(ub(negs)-lb(negs)).*rand(sum(sum(negs)),1);
        end
        if any(pos)
            Island(pos)=lb(pos)+(ub(pos)-lb(pos)).*rand(sum(sum(pos)),1);
        end
    otherwise
        error('flag must be 1, 2 or 3')
end
end

