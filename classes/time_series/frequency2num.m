function fq=frequency2num(fq)
fmap=frequency_map();
if ischar(fq)
    fq=cellstr(fq);
    locs=nan(size(fq));
    for ii=1:numel(fq)
        locs(ii)=find(strcmp(fq{ii},fmap.strings));
    end
    if any(isnan(locs(:)))
        error('some wrong frequency codes found')
    end
    fq=fmap.code(locs);
end
if ~all(ismember(fq,fmap.code))
    error('wrong frequency codes')
end
end