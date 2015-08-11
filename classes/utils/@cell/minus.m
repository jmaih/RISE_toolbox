function c=minus(a,b)
% minus -- overload minus for cell
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

if ischar(a)
    a=cellstr(a);
end

if ischar(b)
    b=cellstr(b);
end

[c,ia]=setdiff(a,b);

[~,isort]=sort(ia);

c=c(isort);

end
