function [a,b]=commute(a,b)
if a.func(1)>b.func(1)
    tmp = a;
    a = b;
    b = tmp;
end
end