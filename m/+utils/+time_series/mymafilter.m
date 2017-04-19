function d=mymafilter(y,q,extend)

if nargin<3
    
    extend=false;
    
end

[T,nvar]=size(y);

N=T-2*q;

wl=2*q+1;

weights=1/(2*q)*ones(wl,1);

weights([1,end])=0.5*weights([1,end]);

weights=weights(:).';

d=zeros(N,nvar);

grab=0:wl-1;

for ii=1:N
    
    grab=grab+1;
    
    d(ii,:)=weights*y(grab,:);
    
end


if extend
        
    d=[d(ones(q,1),:);d;d(N*ones(q,1),:)];
    
    return
    
end

end