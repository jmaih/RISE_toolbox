function c=minus(a,b)
% INTERNAL FUNCTION: overload minus for cell
%

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
