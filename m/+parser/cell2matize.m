function c2m=cell2matize(list,separator)
% INTERNAL FUNCTION: puts list in form (v1|v2|...|vn)
%
% ::
%
%   c2m=cell2matize(list)
%
%   c2m=cell2matize(list,separator)
%
% Args:
%
%    list (cellstr): list of variables
%
%    separator (char | {'|'}): element separating various strings
%
% Returns:
%    :
%
%    - **c2m** [char]: formatted list of variables
%

if nargin<2||isempty(separator)
    
    separator='|';
    
end

c2m=cell2mat(strcat(list(:).',separator));

c2m=['(',c2m(1:end-1),')'];

end