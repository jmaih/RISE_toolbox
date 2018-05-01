function r=resid(obj,varargin)
% resid -- compute the residuals from the steady state
%
% ::
%
%
%   r=resid(obj)
%
%   r=resid(obj,varargin)
%
% Args:
%
%    - **obj** [rise|dsge]: model object
%
%    - **varargin** [pairs of arguments]:
%
% Returns:
%    :
%
%    - **r** [vector|matrix]: residuals
%
% Note:
%
%    - if no output is requested, the residuals are printed on screen
%
% Example:
%
%    See also:

if isempty(obj)
    
    r=cell(0,4);
    
    return
    
end

[~,structural_matrices]=compute_steady_state(obj,varargin{:});

r=full(structural_matrices.user_resids);

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
