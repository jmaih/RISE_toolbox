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
out=regexprep(names,[negative_str,'(\d+)','_(\w+)'],'$2{-$1}');

% variables with leads
positive_str=parser.lead_lag_string(+1);
out=regexprep(out,[positive_str,'(\d+)','_(\w+)'],'$2{+$1}');

end