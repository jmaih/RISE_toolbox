classdef rfvar < svar
    % rfvar Reduced-form VAR modeling
    %
    % rfvar Methods:
    % -----------------
    %
    % check_identification -   H1 line
    % check_optimum -   H1 line
    % draw_parameter -   H1 line
    % estimate -  estimates the parameters of a RISE model
    % forecast -  computes forecasts for rise|dsge|svar|rfvar models
    % get -   H1 line
    % historical_decomposition - Computes historical decompositions of a DSGE model
    % irf -  computes impulse responses for a RISE model
    % isnan -   H1 line
    % load_parameters -   H1 line
    % log_marginal_data_density -   H1 line
    % log_posterior_kernel -   H1 line
    % log_prior_density -   H1 line
    % msvar_priors -   H1 line
    % posterior_marginal_and_prior_densities -   H1 line
    % print_estimation_results -   H1 line
    % prior_plots -   H1 line
    % refresh -  refresh the options of an old object with a newer version of
    % report - assigns the elements of interest to a rise_report.report object
    % rfvar - Reduced-form VAR modeling
    % set -  sets options for RISE models
    % set_solution_to_companion -   H1 line
    % simulate -  simulates a RISE model
    % solve -   H1 line
    % stoch_simul -   H1 line
    % structural_form - finds A structural form given the imposed restrictions
    % template -
    % theoretical_autocorrelations -   H1 line
    % theoretical_autocovariances -   H1 line
    % variance_decomposition -   H1 line
    %
    % rfvar Properties:
    % --------------------
    %
    % identification - information
    % structural_shocks -   information on structural shocks
    % nonlinear_restrictions -   information on nonlinear restrictions
    % constant -   true if VAR has a constant
    % nlags -   number of lags in the VAR
    % legend -   attribute for giving a tag to a specific version of a model
    % endogenous -   information on endogenous variables (names, number, types, etc.)
    % exogenous -   information on exogenous variables (names, number, types, etc.)
    % parameters -   information on parameters (names, number, types, etc.)
    % observables -   information on observable variables (names, number, types, etc.)
    % markov_chains -   information on markov chains, regimes and related items
    % options -   structure holding information on modifiable settings
    % estimation -   information on estimation: posterior maximization and simulation
    % solution -   model solution including steady state, definitions, etc.
    % filtering -   structure holding predicted, updated and smoothed series
    properties(Hidden = true)
        constant_var_data
    end
    properties(SetAccess=protected)
        % identification information
        identification
        % information on structural shocks
        structural_shocks
        % information on nonlinear restrictions
        nonlinear_restrictions
    end
    properties(Access=protected, Hidden = true)
        % the elements below should be for reduced-form vars only
        % hence, the var shall inherit from the structural var
        permute_irf_restrictions=[2,1,3]
        permute_lag_structure_restrictions=[1,2,3]
    end
    
    methods
        function obj=rfvar(varargin)
            obj=obj@svar(varargin{:});
        end
        varargout=structural_form(varargin)
        varargout=check_identification(varargin)
        varargout=solve(varargin)
        varargout=pull_objective(varargin)
    end
    methods(Access=private)
        varargout=translate_restrictions(varargin)
        varargout=set_structural_shocks(varargin)
    end
    methods(Static,Access=private)
        varargout=sort_Q(varargin)
        varargout=var_rotation(varargin)
    end
    methods(Access=protected,Hidden=true)
        varargout=update_estimated_parameter_names(varargin)
        varargout=find_posterior_mode(varargin)
    end
    methods(Static)
        function r=template()
            r=svar.template();
            r.model_class='rfvar';
        end
    end
end

