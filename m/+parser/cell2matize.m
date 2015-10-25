function c2m=cell2matize(list)
% cell2matize -- puts list in form (v1|v2|...|vn)
%
% Syntax
% -------
% ::
%
%   c2m=cell2matize(list)
%
% Inputs
% -------
%
% - **list** [cellstr]: list of variables
%
% Outputs
% --------
%
% - **c2m** [char]: formatted list of variables
%
% More About
% ------------
%
% Examples
% ---------
%
% See also:

c2m=cell2mat(strcat(list,'|'));
c2m=['(',c2m(1:end-1),')'];
end