function warnstate=warnings_fmincon_disable()
% INTERNAL FUNCTION
%

warnstate = warning();% =warning('query','all') %=warning('query') 

warning('off','optim:fmincon:SwitchingToMediumScale')

warning('off','optimlib:fmincon:WillRunDiffAlg')

warning('off','optimlib:fmincon:SwitchingToMediumScaleBecauseNoGrad')

end