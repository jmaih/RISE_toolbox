function [A0Gm_inv,A0Gm]=compute_preconditioner(Gplus01,A0,Q,T)
h=size(Q,1);
n=size(T{1},1);
A0Gm=cell(h,1);
A0Gm_inv=cell(h,1);
for st=1:h
    AQT=A0{st};
    for slead=1:h
        AQT=AQT+Gplus01{st,slead}*T{slead};
    end
    A0Gm{st}=AQT;
    A0Gm_inv{st}=AQT\eye(n);
end