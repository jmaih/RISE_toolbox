function r=weighted_sampling(pop_size,howmany,weights)
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

% samples without repetition
error(nargchk(2,3,nargin))

if nargin<3
    r=utils.optim.uniform_sampling(pop_size,howmany);
else
    if howmany>pop_size
        error([mfilename,':: # requested draws > population size'])
    elseif howmany==pop_size
        r=1:pop_size;
    else
        % first make sure there are at least howmany guys with non-zero
        % weight
        pos=weights>0;
        if sum(pos)<howmany
            error([mfilename,':: not enough positive weights to support sampling'])
        end
        weights=weights/sum(weights);
        cw=cumsum(weights); cw(end)=1;
        pop=1:pop_size;
        r=nan(1,howmany);
        for ii=1:howmany
            index=find(cw>rand,1,'first');
            r(ii)=pop(index);
            pop(index)=[];
            weights(index)=[];
            cw=cumsum(weights); cw(end)=1;
        end
    end
end

