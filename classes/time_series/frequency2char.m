function fq=frequency2char(fq)
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

fmap=frequency_map();
if isnumeric(fq)
    locs=nan(size(fq));
    for ii=1:numel(fq)
        try %#ok<TRYNC>
            locs(ii)=find(fq(ii)==fmap.code);
        end
    end
    if any(isnan(locs(:)))
        error('some wrong frequency codes found')
    end
    fq=char(fmap.strings(locs));
end
if ~all(ismember(fq,fmap.strings))
    error('wrong frequency strings')
end
end
