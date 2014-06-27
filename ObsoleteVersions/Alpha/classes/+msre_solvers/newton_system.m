function [T,G0] = newton_system(T0,Gplus01,A0,Aminus,~,nn,h,~)
% based on the derivatives of the system X+inv(Aplus*X+A0)*Aminus=0. This
% algorithm applies a blind Newton to the equation above.

% compute the G criterion and its derivatives conditional on T0
n2=nn^2;
I2=eye(n2);
Gprime=nan(n2*h);
G0=nan(nn,nn*h);
I1=eye(nn);

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
		Gprime(row,(slead-1)*n2+1:slead*n2)=-kron(LstAm_prime,Lst*Gplus01{st,slead});
	end
    Gprime(row,row)=I2+Gprime(row,row);
end
T0_mat=reshape(T0,nn,nn*h);
T=reshape(T0_mat(:)-Gprime\G0(:),nn,nn,h);
end
