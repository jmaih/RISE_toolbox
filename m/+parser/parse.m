%--- help for inputParser/parse ---
%
% parse(PARSEOBJ,ARGLIST) parses user inputs based on the input scheme.
%    PARSE matches arguments in the comma separated list ARGLIST with
%    input scheme defined by PARSEOBJ.  The inputParser properties Results,
%    UsingDefaults, and Unmatched are set by calling the PARSE methods.  If
%    the parse fails for any reason, the PARSE method throws an error.
% 
%    See also inputParser, inputParser/addOptional, inputParser/addRequired,
%    inputParser/addParameter, inputParser.Results, inputParser.KeepUnmatched, inputParser.UsingDefaults 
%
%    Reference page in Doc Center
%       doc inputParser/parse
%
%