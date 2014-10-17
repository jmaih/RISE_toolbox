function new_name=create_auxiliary_name(name,lead_or_lag)
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

if lead_or_lag
    add_str=parser.lead_lag_string(lead_or_lag);
    new_name=strcat(name,add_str,sprintf('%0.0f',abs(lead_or_lag)));
else
    new_name=name;
end

end