function obj=back_door(obj,attribute,value)
% back_door -- forces the assignment of properties
%
% Syntax
% -------
% ::
%
%   obj=back_door(obj,attribute,value)
%
% Inputs
% -------
%
% - **obj** [rise|dsge|svar|rfvar]: model object
%
% - **attribute** [char|cellstr]: (hierarchical) fields to be updated
%
% - **value** [any]: updating value
%
% Outputs
% --------
%
% - **obj** [generic]: model object
%
% More About
% ------------
%
% Examples
% ---------
%
% See also:

if ischar(attribute)
    
    attribute=cellstr(attribute);
    
end

tmp='obj';

for ii=1:numel(attribute)
    
    tmp=[tmp,'.(attribute{',int2str(ii),'})'];
    
end

tmp=[tmp,'=value;'];

eval(tmp)

% obj.(attribute)=value;

end