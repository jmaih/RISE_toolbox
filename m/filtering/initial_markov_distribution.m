function [PAI00,retcode]=initial_markov_distribution(q_now_lead,ergodic)
if nargin<2
    ergodic=true;
end
nreg=size(q_now_lead,1);
retcode=0;
if ergodic
    A=[eye(nreg)-transpose(q_now_lead);ones(1,nreg)];
    PAI00=A\[zeros(1,nreg),1]';
    if any(isnan(PAI00))||any(PAI00<0)||abs(sum(PAI00)-1)>1e-12
        retcode=308;
    end
else
    PAI00=1/nreg*ones(nreg,1);
end

end
