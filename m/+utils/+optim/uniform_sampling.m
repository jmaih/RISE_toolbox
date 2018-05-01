function r=uniform_sampling(pop_size,howmany)
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

if howmany>pop_size
    error([mfilename,':: # requested draws > population size'])
elseif howmany==pop_size
    r=1:pop_size;
else
    r=1:pop_size;
    for ii=1:howmany
        index=round(ii+rand*(pop_size-ii));    % index=min(np,ceil(ii+rand*(np-ii)));
        swap=r(index);
        r(index)=r(ii);
        r(ii)=swap;
    end
    r=r(1:howmany);
end

