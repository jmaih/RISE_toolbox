function m=symbolic_model(m,var_list)
% H1 line
%
% ::
%
%
% Args:
%
% Returns:
%    :
%
% Note:
%
% Example:
%
%    See also:

if nargin<2
    var_list=parser.input_list();
end
if ischar(var_list)
    var_list=cellstr(var_list);
end

no_word_before='\<'; % no_word_before='(?<!\w+)';

vlist=cell2mat(strcat(var_list(:)','|'));
vlist=vlist(1:end-1);
left_par='\(';
digits='\d+';
optional=',:';
right_par='\)';
pattern=[no_word_before,'(',vlist,')(',left_par,')(',digits,')(',optional,')?(',right_par,')'];

m=regexprep(m,pattern,'$1_$3');