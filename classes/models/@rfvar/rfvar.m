classdef rfvar < svar
    % rfvar class for reduced-form VAR models
    %
    % methods
    % --------
    %
    % - [check_identification](rfvar/check_identification)
    % - [check_optimum](rfvar/check_optimum)
    % - [draw_parameter](rfvar/draw_parameter)
    % - [estimate](rfvar/estimate)
    % - [forecast](rfvar/forecast)
    % - [get](rfvar/get)
    % - [historical_decomposition](rfvar/historical_decomposition)
    % - [irf](rfvar/irf)
    % - [isnan](rfvar/isnan)
    % - [load_parameters](rfvar/load_parameters)
    % - [log_marginal_data_density](rfvar/log_marginal_data_density)
    % - [log_posterior_kernel](rfvar/log_posterior_kernel)
    % - [log_prior_density](rfvar/log_prior_density)
    % - [msvar_priors](rfvar/msvar_priors)
    % - [posterior_marginal_and_prior_densities](rfvar/posterior_marginal_and_prior_densities)
    % - [posterior_simulator](rfvar/posterior_simulator)
    % - [print_estimation_results](rfvar/print_estimation_results)
    % - [prior_plots](rfvar/prior_plots)
    % - [report](rfvar/report)
    % - [rfvar](rfvar/rfvar)
    % - [set](rfvar/set)
    % - [set_solution_to_companion](rfvar/set_solution_to_companion)
    % - [simulate](rfvar/simulate)
    % - [simulation_diagnostics](rfvar/simulation_diagnostics)
    % - [solve](rfvar/solve)
    % - [stoch_simul](rfvar/stoch_simul)
    % - [structural_form](rfvar/structural_form)
    % - [template](rfvar/template)
    % - [theoretical_autocorrelations](rfvar/theoretical_autocorrelations)
    % - [theoretical_autocovariances](rfvar/theoretical_autocovariances)
    % - [variance_decomposition](rfvar/variance_decomposition)
    %
    % properties
    % -----------
    %
    % - [identification] -
    % - [structural_shocks] -
    % - [nonlinear_restrictions] -
    % - [constant] -
    % - [nlags] -
    % - [legend] -
    % - [endogenous] -
    % - [exogenous] -
    % - [parameters] -
    % - [observables] -
    % - [markov_chains] -
    % - [options] -
    % - [estimation] -
    % - [solution] -
    % - [filtering] -
    properties(Hidden = true)
        constant_var_data
    end
    properties(SetAccess=protected)
        identification
        structural_shocks
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
        varargout=update_posterior_simulation_initial_conditions(varargin)
        varargout=initialize_posterior_simulation(varargin)
    end
    methods(Static)
        function r=template()
            r=svar.template();
            r.model_class='rfvar';
            %             r.irf_sign_restrictions=[];
            %             r.irf_zero_restrictions=[];
        end
    end
end

