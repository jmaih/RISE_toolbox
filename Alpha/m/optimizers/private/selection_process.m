function c=selection_process(a,b)

choice=compare_individuals(a,b);
if choice==1
    c=a;
else
    c=b;
end
end