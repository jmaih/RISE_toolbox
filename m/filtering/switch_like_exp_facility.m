function [log_likt,PAI01_tt,retcode]=switch_like_exp_facility(PAI,log_f01,kalman_tol)

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

    retcode=0;
    
    maxf=max(log_f01);
    
    PAI01y=PAI.*exp(log_f01-maxf);
    
    PSTAR=sum(PAI01y);
    
    % Probability updates
    %--------------------
    PAI01_tt=PAI01y/PSTAR;
    
    % Likelihood computation
    %-----------------------
    log_likt=log(PSTAR)+maxf;
    
    if (any(isnan(PAI01_tt))||any(isinf(PAI01_tt))) % log_likt<log(kalman_tol) &&
        
        retcode=306;
               
    end
    
end