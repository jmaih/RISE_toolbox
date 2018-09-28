function o=frwzsp(thetaj,thetaj_bar,pj,iota1,sig)
% thetaj : switching parameter
% thetaj_bar: steady state value
% pj : 1 if the parameter is "partitioned"
% iota1 : 0 if FRWZ perturbation
% sig : perturbation parameter

% % if isa(thetaj_bar,'double') && isnan(thetaj_bar)
% %     % in this way this can be used even in non FRWZ perturbations
% %     %------------------------------------------------------------
% %     thetaj_bar=0;
% %     
% % end

o=(1-(1-iota1)*pj)*thetaj+...
    (1-iota1)*pj*(thetaj_bar+sig*(thetaj-thetaj_bar));

end