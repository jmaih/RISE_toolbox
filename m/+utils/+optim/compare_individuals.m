function c=compare_individuals(a,b)
% INTERNAL FUNCTION
%

c=1;
if (b.violstrength<a.violstrength)||...
        ((b.violstrength==a.violstrength) && b.f<a.f)
    c=2;
end
end