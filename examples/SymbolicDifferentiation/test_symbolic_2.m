clear all
clc
%% create function
func=@(a,b,c,d)[exp(a+2*log(b+c)-a*atan(b*c));exp(-a*atan(b*c));d*exp(a+2*log(b+c))];
%% initialize arguments
a=rise_sad('a');b=rise_sad('b');c=rise_sad('c');d=rise_sad('d');
%% create tree
tree=func(a,b,c,d);
%% compute derivatives... and register the number of calls to each node in the tree 
[dd,references]=diff(tree,{a,b,c,d});
%% now the derivatives can be printed
for irow=1:size(dd,1)
    for jcol=1:size(dd,2)
        disp(['d(',int2str(irow),',',int2str(jcol),')=',char(dd(irow,jcol))])
    end
end