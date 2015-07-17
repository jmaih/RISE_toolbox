function new_name=create_auxiliary_name(name,lead_or_lag,add_prefix)
if nargin<3
    add_prefix=false;
end
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


if ~(ischar(name)||iscellstr(name))
    error('first argument must be a char or a cellstr')
end

if abs(lead_or_lag)>0
    if sign(lead_or_lag)==1
        type='LEAD';
    elseif sign(lead_or_lag)==-1
        type='LAG';
    end
    new_name=sprintf('%s_%0.0f_%s',type,abs(lead_or_lag),name);
else
    if add_prefix
        new_name=sprintf('AUX_%s',name);
    else
        new_name=name;
    end
end

% if lead_or_lag
%     add_str=parser.lead_lag_string(lead_or_lag);
%     new_name=strcat(name,add_str,sprintf('%0.0f',abs(lead_or_lag)));
% else
%     new_name=name;
% end

end