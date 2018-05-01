function m=vectorized_model(m,var_list)
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
vlist=cell2mat(strcat(var_list(:)','|'));
vlist=vlist(1:end-1);

no_word_before='\<';% no_word_before='(?<!\w+)';

left_par='\(';
digits='\d+';
right_par='\)';
underscore='_';
pattern1=[no_word_before,'(',vlist,')(',left_par,')(',digits,')(',right_par,')'];
m=regexprep(m,pattern1,'$1($3,:)');
pattern2=[no_word_before,'(',vlist,')(',underscore,')(',digits,')'];
m=regexprep(m,pattern2,'$1($3,:)');
pattern3='(?<!\.)(/|*|\^)';
m=regexprep(m,pattern3,'.$1');


