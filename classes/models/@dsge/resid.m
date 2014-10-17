function r=resid(obj,varargin)
% H1 line
%
% Syntax
% -------
% ::
%
% Inputs
% -------
%
% Outputs
% --------
%
% More About
% ------------
%
% Examples
% ---------
%
% See also: 

% this function displays the initial residuals of a model (before checking
% that the inital guess actually delivers the steady state)
if isempty(obj)
    r=struct();
    return
end

[obj,ss_and_bgp_start_vals,ssfuncs,retcode]=initial_steady_state(obj,varargin{:});

if retcode
    if obj.options.debug
        utils.error.decipher(retcode)
    end
else
    endo_nbr=obj.endogenous.number(end);
    x=ss_and_bgp_start_vals(1:endo_nbr,:);
    % r is potentially a matrix rather than a vector
    %-----------------------------------------------
    r=ss_residuals(x,ssfuncs.static,ssfuncs.jac_static,obj);
    r=reshape(r,endo_nbr,[]);
    [eqtns_nbr,number_of_regimes]=size(r);
    add_on=repmat('%0.4g ',1,number_of_regimes);
    if nargout==0
        disp(' ')
        for ieqtn=1:eqtns_nbr
            these_resids=num2cell(r(ieqtn,:));
            fprintf(1,['equation #%0.0f : ',add_on,' \n'],ieqtn,these_resids{:});
        end
        clear r
    end
end
