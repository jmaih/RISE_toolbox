function fq=frequency2num(fq)
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

if ischar(fq)
    
    fq=cellstr(fq);
    
end

locs=nan(size(fq));

if iscellstr(fq)
    find_func=@(ii)find(strcmp(fq{ii},fmap.strings));
elseif isnumeric(fq)
    find_func=@(ii)find(fq(ii)==fmap.code);
else
    error('A problem I do not understand')
end

for ii=1:numel(fq)
    
    try %#ok<TRYNC>
        
        locs(ii)=find_func(ii);
        
    end
    
end

if any(isnan(locs(:)))
    
    error('some wrong frequency codes found')
    
end

fq=reshape(fmap.code(locs),size(locs));


if ~all(ismember(fq,fmap.code))
    
    error('wrong frequency codes')
    
end

end