function [T,K0] = functional_iteration(T0,Gplus01,A0,Aminus,~,nn,h,~)

% based on the system X+inv(Aplus*X+A0)*Aminus=0. This
% algorithm reinterprets ??? to derive a New...

% compute the K criterion and its derivatives conditional on T0
K0=nan(nn,nn*h);
T=T0;

for st=1:h % function loop
    AQT=A0{st};
        for slead=1:h
            % AQT=AQT+Q(st,slead)*Aplus{st}*T0(:,:,slead);
            AQT=AQT+Gplus01{st,slead}*T0(:,:,slead);
        end
    K0(:,(st-1)*nn+1:st*nn)=AQT*T0(:,:,st)+Aminus{st};
    T(:,:,st)=-AQT\Aminus{st};
    if any(isnan(T(:)))
        return
    end
end

end