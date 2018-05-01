function A = mergestructures(varargin)
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

% Merge structures with unique fields.
% based on MERGESTRUCT by Loren Shure (Mathworks)

fn = [];
for k = 1:nargin
    try
        fn=[fn;fieldnames(varargin{k})]; %#ok<*AGROW>
    catch MEstruct
        throw(MEstruct)
    end
end

% Make sure the field names are unique.
unique_fields=unique(fn);
isunique=true(1,numel(unique_fields));
for ii=1:numel(unique_fields)
    if sum(strcmp(unique_fields{ii},fn))>1
        isunique(ii)=false;
    end
end
if any(~isunique)
    disp(unique_fields(~isunique))
    error([mfilename,':: Field names must be unique']);
end

% Now concatenate the data from each struct.  Can't use
% structfun since input structs may not be scalar.
c = [];
for k = 1:nargin
    try
        c=[c;struct2cell(varargin{k})];
    catch MEdata
        throw(MEdata);
    end
end

A = orderfields(cell2struct(c,fn,1));
