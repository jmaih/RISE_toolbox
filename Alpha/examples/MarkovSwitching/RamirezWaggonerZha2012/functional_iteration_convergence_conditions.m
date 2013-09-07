function [criterion,aplus,aminus,tau,delta_zero,delta_bk]=functional_iteration_convergence_conditions(model)
aplus=0;
aminus=0;
tau=0;
delta_zero=0;
delta_bk=0;
Aplus=model.Aplus;
Aminus=model.Aminus;
A0=model.A0;
Q=model.Q;
T=model.T;
theta=0;
NumberOfRegimes=model.markov_chains.regimes_number;
for ireg=1:NumberOfRegimes
    aminus=max(aminus,norm(Aminus(:,:,ireg)));
    aplus=max(aplus,norm(Aplus(:,:,ireg)));
    AAT=A0(:,:,ireg);
    for slead=1:NumberOfRegimes
        AAT=AAT+Q(ireg,slead)*Aplus(:,:,slead)*T(:,:,slead);
    end
    tau=max(tau,norm(inv(AAT)));
    theta=max(theta,norm(T(:,:,ireg)));
    delta_zero=max(delta_zero,norm(T(:,:,ireg)-0));
    T_bk=-A0(:,:,ireg)\Aminus(:,:,ireg);
    delta_bk=max(delta_bk,norm(T(:,:,ireg)-T_bk));
end
criterion=tau*aplus*theta./(1-tau*aplus*[delta_zero,delta_bk]);