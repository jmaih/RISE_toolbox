function [T,K0] = functional_iteration_(T0,Aplus,A0,Aminus,Q,nn,h,frwz)
if nargin<8
    frwz=false;
end

% based on the system X+inv(Aplus*X+A0)*Aminus=0. This
% algorithm reinterprets ??? to derive a New...

% compute the G criterion and its derivatives conditional on T0
K0=nan(nn,nn*h);
T=T0;
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
    K0(:,(st-1)*nn+1:st*nn)=AQT*T0(:,:,st)+Aminus(:,:,st);
    T(:,:,st)=-AQT\Aminus(:,:,st);
    if any(any(isnan(T(:,:,st))))
        return
    end
end
end