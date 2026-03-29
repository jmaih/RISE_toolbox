%--- mydeblank.m not found. Showing help for deblank instead. ---
%
% DEBLANK Remove trailing blanks
%    R = DEBLANK(S) removes any trailing whitespace and null characters from
%    the end of S. However, DEBLANK does not remove significant whitespace
%    characters (such as (char(160)).
% 
%    S can be a string array, character array, or a cell array of character
%    vectors. If S is a string array or cell array, then DEBLANK removes 
%    trailing blanks from each element. The output argument R has the same 
%    data type as S.
% 
%    Examples:
%    S = {'MATLAB    ','SIMULINK    '};
%    S = deblank(S)
%    S = 
%       'MATLAB'    'SIMULINK'
% 
%    S = ["Gemini    ","Apollo    ";
%         "ISS       ","Skylab    "];
%    S = deblank(S)
%    S = 
%        "Gemini"    "Apollo"
%        "ISS"       "Skylab"
%        
%    See also ISSPACE, CELLSTR, STRING, STRIP, STRTRIM, PAD.
%
%    Documentation for deblank
%       doc deblank
%
%