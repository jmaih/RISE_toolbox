function g=group(h,varargin)
% INTERNAL FUNCTION
%

% Groups contributions in a structure of ts objects
%
% Args:
%
%    h (struct): structure of ts objects
%    varargin (): same as varargin in ts/group
%
% Returns:
%    :
%
%    - **g** [struct]: structure of ts objects
%
% See also:
%    - ts/group
%

if ~isstruct(h)

    error('this function is written for structures of ts objects')

end

fields=fieldnames(h);

g=struct();

for ifield=1:numel(fields)

    g.(fields{ifield})=group(h.(fields{ifield}),varargin{:});

end

end