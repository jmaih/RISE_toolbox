function [shocks,regimes,PAI,retcode,cfkst,Qt]=forecast(model,y0,ycond,econd,options,regimes)
% FORECAST - conditional forecast for regime-switching models
%
% ::
%
%
%   shocks=forecast(model,y0,ycond,econd,options)
%
%   shocks=forecast(model,y0,ycond,econd,options,regimes)
%
%   [shocks,regimes,retcode,cfkst]=forecast(...)
%
% Args:
%
%    - **model** [struct]: structure with fields
%      - **T** [1 x h cell array]: solution of the model
%      - **sstate** [1 x h cell array]: steady state in each regime
%      - **state_cols** [vector]: location of state variables in the y0 vector
%      below
%      - **Q** [h x h matrix]: initial transition matrix
%      - **Qfunc** [function_handle]: function that takes as input a vector of
%      endogenous variables' values and return a transition matrix
%      - **k** [scalar]: shock horizon beyond the current period
%
%    - **y0** [vector]: initial conditions
%
%    - **ycond** [struct]: structure with fields
%      - **data** [3-dimensional array]: The first page is the mean, the
%      second is the lower bound, the third is the upper bound
%      - **pos** [empty|scalar|vector]: location of the conditioning variables
%      in the state vector
%
%    - **econd** []: structure with fields
%      - **data** [3-dimensional array]: The first page is the mean, the
%      second is the lower bound, the third is the upper bound
%      - **pos** [empty|scalar|vector]: location of the conditioning shocks
%      in the state vector
%
%    - **options** [struct]: structure with options of interest. See details
%    in utils.forecast.rscond.density_shocks and utils.forecast.rscond.form_system
%
%    - **regimes** [empty|vector]: regimes for each step
%
% Returns:
%    :
%
%    - **shocks** [empty|2 or 3-dimensional array]: nshocks x horizon x ndraws
%
%    - **regimes** [vector]: regimes visited over the forecast horizon
%
%    - **PAI** [h x nsteps matrix]: evolution of choice probabilities
%
%    - **retcode** [scalar]: 0 if no problem
%
%    - **cfkst** [n x (nsteps+1) matrix]: conditional forecasts
%
%    - **Qt** [h x h x nsteps array]: time series of transition matrices
%
% Note:
%
% Example:
%
%    See also: RSCF.LOOP_FORECAST

if nargin==0
    
    [s1,def1]=utils.forecast.rscond.form_system();
    
    [s2,def2]=utils.forecast.rscond.density_shocks();
    
    remove_duplicates();
    
    if nargout
        
        shocks=[def1;def2];
        
    else
        
        disp_defaults([def1;def2]);
        
    end
    
    return
    
else
    narginchk(5,6)
    
    if nargin<6
        
        regimes=[];
        
    end
    
    nsteps=options.nsteps;
    
    [M,~,regimes,PAI,~,Qt]=utils.forecast.rscond.stochastic_impact_cumulator(model,y0,nsteps,ycond.pos,...
        econd.pos,regimes);
        
    % expand the data and chop them if necessary
    %--------------------------------------------
    [R_S,mu_lb_ub]=utils.forecast.rscond.form_system(M,ycond,econd,nsteps,...
        utils.miscellaneous.reselect_options(options,'utils.forecast.rscond.form_system'));
    
    % find shocks and reshape them
    %------------------------------
    [shocks,retcode]=utils.forecast.rscond.density_shocks(R_S,mu_lb_ub,model.nshocks,[],...
        utils.miscellaneous.reselect_options(options,'utils.forecast.rscond.density_shocks'));
    
    % Compute forecasts
    %--------------------
    cfkst=[];
    
    if ~retcode && nargout>3
        
        cfkst=utils.forecast.rscond.condition_on_shocks(model.T,model.sstate,y0,regimes,...
            model.state_cols,model.k,options.nsteps,shocks,options.debug);
        
        is_hard=~isempty(ycond.data) && max(max(abs(ycond.data(:,:,3)-ycond.data(:,:,2))))<1e-8;
        
        if is_hard && options.debug
            
            cy=size(ycond.data,2);
            
            max_discrep=max(max(abs(cfkst(ycond.pos,1+(1:cy))-ycond.data(:,:,1))));
            
            fprintf('Checking hard conditions max discrep = %0.8f\n\n',max_discrep)
        
        end
        
    end
    
end

    function remove_duplicates()
        
        fields1=fieldnames(s1);
        
        fields2=fieldnames(s2);
        
        bad=ismember(fields1,fields2);
        
        bad_list=fields1(bad);
        
        s1=rmfield(s1,bad_list);
        
        def1=def1(~bad,:);
        
    end

end