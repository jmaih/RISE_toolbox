function obj=setup_general_restrictions(obj)
% setup_general_restrictions - sets up the general nonlinear restrictions
%
% ::
%
%
%	obj=setup_general_restrictions(obj)
%
% Args:
%
%    - **obj** [rise|dsge|rfvar|svar]: RISE object
%
% Returns:
%    :
%
%    - **obj** [rise|dsge|rfvar|svar]: RISE object
%
% Note:
%
%    - If present, the general restrictions should have been passed at an earlier stage
%   	See help for generic/estimate for further information
%
% Example:
%
%    See also:

% this can be useful also for checking restrictions such as explosive VARs,
% etc. and not just restrictions on DSGEs
nobj=numel(obj);

for ii=1:nobj
    
    if ~isempty(obj(ii).options.estim_general_restrictions)

		[genrestr,vargs]=utils.code.user_function_to_rise_function(...
                    obj(ii).options.estim_general_restrictions);
        
        obj(ii).general_restrictions_data=@(x)genrestr(x,vargs{:});
            
    end
    
    % check the filtering level required under estimation
    %-----------------------------------------------------
    if ~isempty(obj(ii).general_restrictions_data)
        
        kflev=load_kalman_filtering_level('estim_general_restrictions');
                
        obj(ii).options.kf_filtering_level=kflev;
        
    end
    
    % all the priors have already been set. Check whether there are
    % endogenous priors
    if ~isempty(obj.options.estim_endogenous_priors)
        
        kflev=load_kalman_filtering_level('estim_endogenous_priors');
        
        obj(ii).options.kf_filtering_level=max(kflev,...
            obj(ii).options.kf_filtering_level);
        
    end
    
end

    function kflev=load_kalman_filtering_level(field)
        
        genrest_basic=obj(ii).options.(field);
        
        if iscell(genrest_basic)
            
            genrest_basic=genrest_basic{1};
            
        end
        
        req=genrest_basic();
        
        kflev=req.kf_filtering_level;
        
    end

end