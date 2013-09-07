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
if ~frwz
    % Create a variable that will hold Aplus*T and delete Aplus
    AplusT0=Aplus;
end
% Get the fwl and trim Aplus once and for all
% Get the bkl and trim Aminus once and for all
fwl=any(Aplus(:,:,1));
bkl=any(Aminus(:,:,1));
for istate=2:h
    fwl=fwl|any(Aplus(:,:,istate));
    bkl=bkl|any(Aminus(:,:,istate));
end
Aplus=Aplus(:,fwl,:);
Aminus=Aminus(:,bkl,:);
for st=1:h % function loop
    if frwz
        AQT=0;
        for slead=1:h
            AQT=AQT+Q(st,slead)*T0(:,:,slead);
        end
        AQT=A_times_X(Aplus(:,:,st),AQT,fwl)+A0(:,:,st);
    else
        AQT=A0(:,:,st);
        for slead=1:h
            if st==1
                % this is so that those computations are done only once
                AplusT0(:,:,slead)=A_times_X(Aplus(:,:,slead),T0(:,:,slead),fwl);
            end
            AQT=AQT+Q(st,slead)*AplusT0(:,:,slead);
        end
        % there is nothing to clear in this case since we sill need Aplus
        % and that we need AplusT0 throughout the loop
    end
    Lst=AQT\I1;
    if any(any(isnan(Lst)))
        T=nan(size(T0));
        return
    end
    row=(st-1)*n2+1:st*n2;
    LstAm_prime=transpose(...
        X_times_A(Lst,Aminus(:,:,st),bkl)...
        );
    G0(:,(st-1)*nn+1:st*nn)=T0(:,:,st)+transpose(LstAm_prime);
    if frwz
        CC(row,:)=kron(Q(st,:),kron(LstAm_prime,...
            X_times_A(Lst,Aplus(:,:,st),fwl)...
            ));
    else
        for jj=1:h
            col=(jj-1)*n2+1:jj*n2;
            CC(row,col)=Q(st,jj)*kron(LstAm_prime,...
                X_times_A(Lst,Aplus(:,:,jj),fwl)...
                );
        end
    end
    CC(row,row)=CC(row,row)-I2;
end
DELTA=reshape(CC\G0(:),nn,nn,h);
T=T0+DELTA;

end
