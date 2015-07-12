function p=derivative_location(a,p)

if nargin<2
    p=0;
end
o=numel(a);

if o==1
    p=p+a;
else
    a_1=a(1)-1;
    if a_1
        if a(1)==a(end)
            p0=utils.gridfuncs.my_nchoosek(a(1)+o-1,o);
        else
            p0=utils.gridfuncs.my_nchoosek(a_1+o-1,o);
            p0=utils.gridfuncs.derivative_location(a(2:end),p0);
        end
    else
        p0=1;
    end
    p=p+p0;
end

end