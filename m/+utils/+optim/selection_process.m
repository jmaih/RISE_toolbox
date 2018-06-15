function c=selection_process(a,b)
% INTERNAL FUNCTION
%

choice=utils.optim.compare_individuals(a,b);
if choice==1
    c=a;
else
    c=b;
end
end