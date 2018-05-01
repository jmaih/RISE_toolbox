function ind=stationary_index(obj,too_small)

% STATIONARY_INDEX - index for stationary variables
%
% ::
%
%
%   ind=STATIONARY_INDEX(obj)
%
%   ind=STATIONARY_INDEX(obj,too_small)
%
% Args:
%
%    - **obj** [rise|dsge]: model object
%
%    - **too_small** [numeric|{1e-10}]: cutoff criterion
%
% Returns:
%    :
%
%    - **ind** [logical]: true if variable is stationary
%
% Note:
%
% Example:
%
%    See also:

if isempty(obj)
    
    ind=struct();
    
    return
    
end

if nargin<2
    
    too_small=1e-10;
    
end

% solution is calculated (hence the bgp) before transforming it for
% log_expanded variables and so, in this case, their bgp is computed in the
% same way as for level variables. Hence log_expanded variables do not belong
% in the list below

if isempty(obj.solution)||~isfield(obj.solution,'bgp')
    
    error('The model needs to be solved first')
    
end

n=obj.endogenous.number;

ind=true(1,n);

if ~obj.options.solve_bgp
    
    return
    
end

bgp=cell2mat(obj.solution.bgp);

checklog=@(x)~any(abs(bgp(x,:)-1)>too_small);

checklev=@(x)~any(abs(bgp(x,:))>too_small);

log_vars=obj.endogenous.is_log_var;

for ivar=1:n
    
    ind(ivar)=if_then_else(log_vars(ivar),checklog(ivar),checklev(ivar));
    
end

end