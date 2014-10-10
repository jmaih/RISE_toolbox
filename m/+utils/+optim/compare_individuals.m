function c=compare_individuals(a,b)
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

c=1;
if (b.violstrength<a.violstrength)||...
        ((b.violstrength==a.violstrength) && b.f<a.f)
    c=2;
end
end