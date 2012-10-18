function [T,K0] = functional_iteration(T0,Aplus,A0,Aminus,Q,nn,h,frwz)
if nargin<8
    frwz=false;
end

% based on the system X+inv(Aplus*X+A0)*Aminus=0. This
% algorithm reinterprets ??? to derive a New...

fwl=any(Aplus(:,:,1));
for istate=2:h
    fwl=fwl|any(Aplus(:,:,istate));
end

% compute the G criterion and its derivatives conditional on T0
K0=nan(nn,nn*h);
T=T0;
if ~frwz
    % Create a variable that will hold Aplus*T and delete Aplus
    AplusT0=Aplus;
end
% now we trim Aplus to save on memory
Aplus=Aplus(:,fwl,:);
for st=1:h % function loop
    if frwz
        AQT=0;
        for slead=1:h
            AQT=AQT+Q(st,slead)*T0(:,:,slead);
        end
        AQT=A_times_X(Aplus(:,:,st),AQT,fwl)+A0(:,:,st); % NB: Aplus us trimmed
    else
        AQT=A0(:,:,st);
        for slead=1:h
            if st==1
                % this is so that those computations are done only once
                AplusT0(:,:,slead)=A_times_X(Aplus(:,:,slead),T0(:,:,slead),fwl); % NB: Aplus us trimmed
            end
            AQT=AQT+Q(st,slead)*AplusT0(:,:,slead);
        end
        clear Aplus
    end
    K0(:,(st-1)*nn+1:st*nn)=AQT*T0(:,:,st)+Aminus(:,:,st);
    T(:,:,st)=-AQT\Aminus(:,:,st);
    if any(any(isnan(T(:,:,st))))
        return
    end
end
end