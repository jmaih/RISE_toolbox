function fq=frequency2char(fq)
fmap=frequency_map();
if isnumeric(fq)
    locs=nan(size(fq));
    for ii=1:numel(fq)
        locs(ii)=find(fq(ii)==fmap.code);
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
