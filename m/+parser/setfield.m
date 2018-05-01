function S = setfield(S,name,value)
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

%UNTITLED10 Summary of this function goes here
%   Detailed explanation goes here

siz=size(S);

S=S(:);

fnames=fieldnames(S);

S=struct2cell(S);

ncols=size(S,2);

if ischar(value)
    value=cellstr(value);
elseif isa(value,'double')||isa(value,'logical')
    value=num2cell(value);
end
value=value(:)';

if numel(value)==ncols
elseif numel(value)==1 && ncols >1
    value=repmat(value,1,ncols);
else
    error('mismatch between the number of cases of S and the number of values to be added to the structure')
end

S=[S;value];

fnames=[fnames(:);name];

S=cell2struct(S,fnames,1);

S=reshape(S,siz);

end

