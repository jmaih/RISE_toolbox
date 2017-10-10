function mavg=moving_average(B,nz,order)

T=vartools.companion(B,[],nz);

mavg=cell(1,order);

ncomp=size(T,1);

n=size(B,1);

Tpower=eye(ncomp);

for io=1:order+1
    
    mavg{io}=Tpower(1:n,1:n);
    
    Tpower=T*Tpower;
    
end

end