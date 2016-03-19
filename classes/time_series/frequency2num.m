function fq=frequency2num(fq)
% H1 line
%
% Syntax
% -------
% ::
%
% Inputs
% -------
%
% Outputs
% --------
%
% More About
% ------------
%
% Examples
% ---------
%
% See also: 

fmap=frequency_map();
if ischar(fq)
    fq=cellstr(fq);
    locs=nan(size(fq));
    for ii=1:numel(fq)
        try %#ok<TRYNC>
        locs(ii)=find(strcmp(fq{ii},fmap.strings));
        end
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