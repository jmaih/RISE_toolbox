function c2m=cell2matize(list)
% cell2matize -- puts list in form (v1|v2|...|vn)
%
% ::
%
%
%   c2m=cell2matize(list)
%
% Args:
%
%    - **list** [cellstr]: list of variables
%
% Returns:
%    :
%
%    - **c2m** [char]: formatted list of variables
%
% Note:
%
% Example:
%
%    See also:

c2m=cell2mat(strcat(list,'|'));
c2m=['(',c2m(1:end-1),')'];
end