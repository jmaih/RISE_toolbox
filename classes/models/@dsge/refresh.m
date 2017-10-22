function varargout=refresh(obj)
% REFRESH - refresh the options of an old object with a newer version of
%   the software
%
% More About
% ------------
%
% - REFRESH is the same as RISE_GENERIC.REFRESH except that it also
% refreshes parts of the system that are specific to DSGE or RISE objects.
%
% Examples
% ---------
%
% See also: RISE_GENERIC.REFRESH

if ~isempty(obj)
    
    nregs=obj.markov_chains.regimes_number;
    
    if nregs>1
        
        if size(obj.exogenous.shock_horizon,1)==1
            
            obj.exogenous.shock_horizon=obj.exogenous.shock_horizon(ones(nregs,1),:);
            
        end
        
    end
    
end

[varargout{1:nargout}]=refresh@generic(obj);

end