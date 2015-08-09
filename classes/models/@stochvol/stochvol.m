classdef  stochvol < rfvar
    % stochvol class for time-varying parameter models
    %
    % methods
    % --------
    %
    % - [check_identification](stochvol/check_identification)
    % - [check_optimum](stochvol/check_optimum)
    % - [draw_parameter](stochvol/draw_parameter)
    % - [estimate](stochvol/estimate)
    % - [forecast](stochvol/forecast)
    % - [get](stochvol/get)
    % - [historical_decomposition](stochvol/historical_decomposition)
    % - [irf](stochvol/irf)
    % - [isnan](stochvol/isnan)
    % - [load_parameters](stochvol/load_parameters)
    % - [log_marginal_data_density](stochvol/log_marginal_data_density)
    % - [log_posterior_kernel](stochvol/log_posterior_kernel)
    % - [log_prior_density](stochvol/log_prior_density)
    % - [msvar_priors](stochvol/msvar_priors)
    % - [posterior_marginal_and_prior_densities](stochvol/posterior_marginal_and_prior_densities)
    % - [print_estimation_results](stochvol/print_estimation_results)
    % - [prior_plots](stochvol/prior_plots)
    % - [report](stochvol/report)
    % - [set](stochvol/set)
    % - [set_solution_to_companion](stochvol/set_solution_to_companion)
    % - [simulate](stochvol/simulate)
    % - [solve](stochvol/solve)
    % - [stoch_simul](stochvol/stoch_simul)
    % - [stochvol](stochvol/stochvol)
    % - [structural_form](stochvol/structural_form)
    % - [template](stochvol/template)
    % - [theoretical_autocorrelations](stochvol/theoretical_autocorrelations)
    % - [theoretical_autocovariances](stochvol/theoretical_autocovariances)
    % - [variance_decomposition](stochvol/variance_decomposition)
    %
    % properties
    % -----------
    %
    % - [time_varying_parameters] -
    % - [random_walk_parameters] -
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
    properties(Hidden=true)
        endogenous_positions
        exogenous_positions
        parameters_positions
    end
    properties(SetAccess=protected)
        time_varying_parameters
        random_walk_parameters
    end
    methods
        function obj=stochvol(varargin)
            obj=obj@rfvar(varargin{:});
            if ~isempty(varargin)
                obj.time_varying_parameters=logical(varargin{1}.time_varying_parameters);
                obj.random_walk_parameters=logical(varargin{1}.random_walk_parameters);
                if ~any(obj.time_varying_parameters)
                    error('no time varying parameters and no stochastic volatility, please use the more efficient "rfvar" object instead')
                end
                % redefine endogenous,exogenous,parameters, etc.
                %-----------------------------------------------
                obj=stochvol.redo_declarations(obj);
                % link the parameters to the structural matrices
                %-----------------------------------------------
                [obj.param_to_mat_links,obj.parameter_values]=...
                    vartools.parameters_to_matrices(...
                    obj.param_template,obj.parameters.name,...
                    obj.markov_chains.regimes_number);
            end
        end
    end
    methods(Sealed)
    end
    methods(Sealed,Hidden=true)
        varargout=simulation_engine(varargin)
    end
    methods(Static,Hidden=true)
        varargout=redo_declarations(varargin)
        varargout=locate_variables_blocks(varargin)
        varargout=form_parameter_matrices(varargin)
        varargout=format_blocks(varargin)
        % % % % % % % % % % % % % %         varargout=gibbs_sampler(varargin)
    end
    methods(Access=protected)
        varargout=setup_linear_restrictions(varargin)
    end
    methods(Static)
        function r=template()
            r=rfvar.template();
            r.model_class='stochvol';
            % lags,volatility,correlations
            r.time_varying_parameters=[true,true,false];
            r.random_walk_parameters=true(1,3);
        end
    end
end

