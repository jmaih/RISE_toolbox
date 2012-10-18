function [ishigami,bounds]=satelli_sobol95()

ishigami=@gfunc;

bounds=repmat([0,1],8,1);

end

function f=gfunc(x)
debug=false;
a=[4.5,4.5,1,0,1,9,0,9];
f=1;
for ii=1:numel(a)
    gi=(abs(4*x(ii,:)-2)+a(ii))/(1+a(ii));
    f=f.*gi;
end
if debug
    f=f(ones(10,1),:);
end
end