function x=reformat_restriction(x,convertfunc1,convertfunc2) %#ok<INUSD>
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

% the first conversion function finds the equation
% the second conversion function finds the column, which may be
% either an endogenous variable or a shock
if nargin<3
    convertfunc2=convertfunc1; %#ok<NASGU>
end
cellmode=iscellstr(x);
if ~cellmode
    x=cellstr(x);
end
nrest=numel(x);
for ix=1:nrest
    x{ix}(isspace(x{ix}))=[];
end
% add date 0 to undated variables i.e. e.g. Y@R-->Y@R{0}
%-------------------------------------------------------
x=regexprep(x,'(\w+)(@\w+)(?!({|\w+))','$1$2{0}');
% remove minus or plus signs in curly braces
%-------------------------------------------
x=regexprep(x,'({)(\-|\+)(\d+})','$1$3');
% put in a standard form a2(v1,v2)
%---------------------------------
x=regexprep(x,'(\w+)@(\w+){(\d+)}','a$3($1,$2)');
% put in the form alag_eqtn_var from either alag(1,v) or
% alag(v1,v2)
%-------------------------------------------------------
x=regexprep(x,'(a\w+)\((\d+),(\w+)\)','$1_$2_${convertfunc2($3)}');
x=regexprep(x,'(a\w+)\((\w+),(\w+)\)','$1_${convertfunc1($2)}_${convertfunc2($3)}');
if ~cellmode
    x=char(x);
end
end
