function [T,G0] = newton_system(T0,Aplus,A0,Aminus,Q,nn,h,frwz)
% based on the derivatives of the system X+inv(Aplus*X+A0)*Aminus=0. This
% algorithm applies a blind Newton to the equation above.
if nargin<8
    frwz=false;
end

% compute the G criterion and its derivatives conditional on T0
n2=nn^2;
In2=eye(n2);
Gprime=nan(n2*h);
G0=nan(nn,nn*h);
In=eye(nn);
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
    f_id=(st-1)*n2+1:st*n2;
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
    AAQTi=AQT\In;
    if any(any(isnan(AAQTi)))
        T=nan(size(T0));
        return
    end
    AAQTiAm_prime=transpose(...
        X_times_A(AAQTi,Aminus(:,:,st),bkl));
    G0(:,(st-1)*nn+1:st*nn)=T0(:,:,st)+transpose(AAQTiAm_prime);
    if frwz
        Gprime(f_id,:)=-kron(Q(st,:),kron(AAQTiAm_prime,...
            X_times_A(AAQTi,Aplus(:,:,st),fwl)...
            ));
    else
        for jj=1:h
            col=(jj-1)*n2+1:jj*n2;
            Gprime(f_id,col)=-Q(st,jj)*kron(AAQTiAm_prime,...
                X_times_A(AAQTi,Aplus(:,:,jj),fwl)...
                );
        end
    end
    Gprime(f_id,f_id)=In2+Gprime(f_id,f_id);
end
T0_mat=reshape(T0,nn,nn*h);
T=reshape(T0_mat(:)-Gprime\G0(:),nn,nn,h);
end
