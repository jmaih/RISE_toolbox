function [T,G0] = newton_kronecker(T0,Aplus,A0,Aminus,Q,nn,h,frwz)
% alternative newton algorithm based on the expansion of the system
% X+inv(Aplus*X+A0)*Aminus=0.
if nargin<8
    frwz=false;
end

% compute the G criterion
n2=nn^2;
I1=eye(nn);
I2=eye(n2);
G0=nan(nn,nn*h);
CC=zeros(n2*h);
for st=1:h % function loop
    if frwz
        AQT=0;
        for slead=1:h
            AQT=AQT+Q(st,slead)*T0(:,:,slead);
        end
        AQT=Aplus(:,:,st)*AQT+A0(:,:,st);
    else
        AQT=A0(:,:,st);
        for slead=1:h
            AQT=AQT+Q(st,slead)*Aplus(:,:,slead)*T0(:,:,slead);
        end
    end
    Lst=AQT\I1;
    if any(any(isnan(Lst)))
        T=nan(size(T0));
        return
    end
    row=(st-1)*n2+1:st*n2;
    LstAm_prime=transpose(Lst*Aminus(:,:,st));
    G0(:,(st-1)*nn+1:st*nn)=T0(:,:,st)+transpose(LstAm_prime);
    if frwz
        CC(row,:)=kron(Q(st,:),kron(LstAm_prime,Lst*Aplus(:,:,st)));
    else
        for jj=1:h
            col=(jj-1)*n2+1:jj*n2;
            CC(row,col)=Q(st,jj)*kron(LstAm_prime,Lst*Aplus(:,:,jj));
        end
    end
    CC(row,row)=CC(row,row)-I2;
end
DELTA=reshape(CC\G0(:),nn,nn,h);
T=T0+DELTA;

end
