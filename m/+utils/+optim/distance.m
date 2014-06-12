function ed=distance(a,b)
if isempty(a)||isempty(b)
    ed=inf;
else
    dev=bsxfun(@minus,a,b);
    ed=sqrt(sum(dev.^2,1));
end
