function [PAI00,retcode]=initial_markov_distribution(q_now_lead,flag)
nreg=size(q_now_lead,1);
retcode=0;
switch flag
    case{0,'ergodic'}
        A=[eye(nreg)-transpose(q_now_lead);ones(1,nreg)];
        PAI00=A\[zeros(1,nreg),1]';
        if any(isnan(PAI00))||any(PAI00<0)||abs(sum(PAI00)-1)>1e-12
            retcode=308;
        end
    case{1,'diffuse'}
        PAI00=1/nreg*ones(nreg,1);
    otherwise
        error([mfilename,':: unkown '])
end

