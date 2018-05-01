function obj=back_door(obj,attribute,value)
% back_door -- forces the assignment of properties
%
% ::
%
%
%   obj=back_door(obj,attribute,value)
%
% Args:
%
%    - **obj** [rise|dsge|svar|rfvar]: model object
%
%    - **attribute** [char|cellstr]: (hierarchical) fields to be updated
%
%    - **value** [any]: updating value
%
% Returns:
%    :
%
%    - **obj** [generic]: model object
%
% Note:
%
% Example:
%
%    See also:

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