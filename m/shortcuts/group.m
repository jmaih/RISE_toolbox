function g=group(h,varargin)
% GROUP -- groups contributions in a structure of ts objects
%
% ::
%
%
% Args:
%
%    - **h** [struct]: structure of ts objects
%
%    - **varargin** []: same as varargin in TS/GROUP
%
% Returns:
%    :
%
%    - **g** [struct]: structure of ts objects
%
% Note:
%
% Example:
%
%    See also: TS/GROUP

if ~isstruct(h)
    
    error('this function is written for structures of ts objects')
    
end

fields=fieldnames(h);

g=struct();

for ifield=1:numel(fields)
    
    g.(fields{ifield})=group(h.(fields{ifield}),varargin{:});
    
end

end