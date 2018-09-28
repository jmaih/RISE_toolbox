function warnstate=warnings_disable()
% INTERNAL FUNCTION
%

warnstate = warning();% =warning('query','all') %=warning('query') 

warning('off','MATLAB:nearlySingularMatrix')

warning('off','MATLAB:illConditionedMatrix')

warning('off','MATLAB:singularMatrix')

end