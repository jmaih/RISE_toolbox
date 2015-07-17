function obj=add_to_routines(obj,varargin)
% add_to_routines - adds routines items
%
% Syntax
% -------
% ::
%
%   obj=add_to_routines(obj,item1_name,item1_val,...)
%
% Inputs
% -------
%
% - **obj** [rise|dsge]: model object
%
% - **varargin** : arguments coming in pairs
%
% Outputs
% --------
%
% - **obj** [rise|dsge]: model object
%
% More About
% ------------
%
% Examples
% ---------
%
% See also:

n=length(varargin);
if rem(n,2)
    error('arguments must come in pairs')
end

names=varargin(1:2:end);

vals=varargin(2:2:end);

if isempty(obj.routines)
    obj.routines=struct();
    obj.online_routines=struct();
end
for item=1:numel(names)
    name=names{item};
    if ~ischar(name)
        error('in pairwise inputs the first is always a char')
    end
    obj.routines.(names{item})=vals{item};
    obj.online_routines.(names{item})=vals{item};
end
end