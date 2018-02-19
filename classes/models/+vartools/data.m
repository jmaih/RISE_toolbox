function kdata = data(mydata,nlags,constant)

[nvars,Tr] = size(mydata);

row_wise=true;
%
X = another_embed(mydata,nlags+1,row_wise);

Y = X(1:nvars,:); % mydata(:,nlags+1:end);

X = X(nvars+1:end,:);

if constant
    
    X=[ones(1,Tr-nlags);X];
    
end
%
K = size(X,1);

T = Tr - nlags;

kdata=struct('Y',Y,'X',X,'nvars',nvars,'K',K,'T',T,...
    'constant',constant,'nlags',nlags);

end

function b=another_embed(a,n,row_wise)

if nargin<3
    
    row_wise=false;
    
end

if row_wise
    
    a=a.';
    
end

b=cell(1,n);

for icol=1:n
    
    b{icol}=a((n:end)-icol+1,:);
    
end

b=cell2mat(b);

if row_wise
    
    b=b.';
    
end

end
