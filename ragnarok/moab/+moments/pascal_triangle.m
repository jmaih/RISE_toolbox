function pt=pascal_triangle(order)

pt=cell(order+1,1);

pt{1}=1;

for ii=2:order+1
    
    pt{ii}=main_engine(pt{ii-1});
    
end

end

function next=main_engine(prev)

n=numel(prev);

next=zeros(1,n+1);

next([1,n+1])=1;

for ii=2:n
    
    next(ii)=prev(ii-1)+prev(ii);
    
end

end