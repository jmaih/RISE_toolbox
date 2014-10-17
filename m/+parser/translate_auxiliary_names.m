function out=translate_auxiliary_names(names)
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


% variables with lags
negative_str=parser.lead_lag_string(-1);
out=regexprep(names,['(\w+)',negative_str,'(\d+)'],'$1{-$2}');

% variables with leads
positive_str=parser.lead_lag_string(+1);
out=regexprep(out,['(\w+)',positive_str,'(\d+)'],'$1{+$2}');

end