function c=selection_process(a,b)
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
% Description
% ------------
%
% Examples
% ---------
%
% See also: 

choice=utils.optim.compare_individuals(a,b);
if choice==1
    c=a;
else
    c=b;
end
end