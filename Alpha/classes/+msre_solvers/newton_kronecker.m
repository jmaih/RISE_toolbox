function [T,G0] = newton_kronecker(T0,Gplus01,A0,Aminus,~,nn,h,~)
%function [T,G0] = newton_kronecker(T0,Gplus01,A0,Aminus,Q,nn,h,frwz)
% alternative newton algorithm based on the expansion of the system
% X+inv(Aplus*X+A0)*Aminus=0.

% compute the G criterion
n2=nn^2;
I1=eye(nn);
I2=eye(n2);
G0=nan(nn,nn*h);
CC=zeros(n2*h);

for st=1:h % function loop
    AQT=A0{st};
    for slead=1:h
        % AQT=AQT+Q(st,slead)*Aplus{st}*T0(:,:,slead);
        AQT=AQT+Gplus01{st,slead}*T0(:,:,slead);
    end
    Lst=AQT\I1;
    if any(isnan(Lst(:)))
        T=nan(size(T0));
        return
    end
    row=(st-1)*n2+1:st*n2;
    LstAm=Lst*Aminus{st};
    G0(:,(st-1)*nn+1:st*nn)=T0(:,:,st)+LstAm;
	LstAm_prime=transpose(LstAm); clear LstAm 
	for slead=1:h
		CC(row,(slead-1)*n2+1:slead*n2)=kron(LstAm_prime,Lst*Gplus01{st,slead});
	end
    CC(row,row)=CC(row,row)-I2;
end
DELTA=reshape(CC\G0(:),nn,nn,h);
T=T0+DELTA;

end
