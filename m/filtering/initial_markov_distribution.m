function [PAI00,retcode]=initial_markov_distribution(q_now_lead,ergodic)
% H1 line
%
% ::
%
%
% Args:
%
% Returns:
%    :
%
% Note:
%
% Example:
%
%    See also:

if nargin<2
    
    ergodic=true;
    
end

nreg=size(q_now_lead,1);

retcode=0;

Qdiag=diag(diag(q_now_lead));

if max(abs(Qdiag(:)-q_now_lead(:)))<1e-12||~ergodic
    
    PAI00=1/nreg*ones(nreg,1);
    
else
    
    PAI00=[eye(nreg)-q_now_lead'
        ones(1,nreg)]\[zeros(nreg,1)
        1];
    
    % safeguard against knife-edge cases e.g. PAI(1)=1 and PAI(2)=1e-17
    %------------------------------------------------------------------
    very_small=abs(PAI00)<1e-9;
    
    if any(very_small)
        
        PAI00(very_small)=0;
        
        PAI00=PAI00/sum(PAI00);
        
    end
    
    if any(isnan(PAI00))||any(PAI00<0)||abs(sum(PAI00)-1)>1e-12
        
        retcode=308;
        
    end
    
end


end
