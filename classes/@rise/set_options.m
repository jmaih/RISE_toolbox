% function obj=set_options(obj,varargin)
%
% for ii=1:numel(obj)
%     [obj(ii).options,missing]=mysetfield(obj(ii).options,varargin{:});
%     if ~isempty(missing)
%         disp([mfilename,':: the following variables are missing'])
%         disp(missing)
%     end
% end

function obj=set_options(obj,varargin)

init=isempty(obj.options);

if init
    % create a dummy object to read the options from various
    % objects, functions and sub-functions
    dum_obj=rise.empty(0,1);
    method_list=methods(dum_obj);
    restricted_list={'Contents','rise','set_options','minimax',...
        'check_optimum','posterior_marginal_and_prior_densities',...
        'print_estimates','print_estimation_results','print_solution','prior_plots',...
        'simulation_diagnostics','theoretical_autocovariances'};
    method_list(ismember(method_list,restricted_list))=[];
    % add some private methods to the list
    method_list=[method_list;'load_functions'];

    myoptions=struct();
    for ii=1:numel(method_list)
        %                 if  nargout(['rise>rise.',method_list{ii}])~=0
        myoptions=mergestructures(myoptions,...
            dum_obj.(method_list{ii})...
            );
        %                 end
    end
    
    % this is needed later to load the functions
    refresh_functions=true;
    
    % Create the main options
    update_random_number_stream=true;
    expansion_order=1;
    % shocks are allowed to be anticipated.
    shock_properties=struct('name',{obj.varexo.name}',...
        'StandardDeviation',num2cell(nan(obj.NumberOfExogenous,1)),...
        'horizon',num2cell(expansion_order*ones(obj.NumberOfExogenous,1)));
    rows=4;cols=3;
    MainOptions={
         % set output folder name
		'results_folder',obj.filename %  obj.output_folder_name=obj.filename;
        % ====== function handle that solves the steady state ======
        'use_steady_state_model',true
        % ====== forecast & time series options ======
        'data',rise_time_series
        'data_demean',false
        'real_time',true
        % if real_time is set to true, the horizon of the conditional
        % information is interpreted as real time data. When real_time is
        % false, the horizon of the conditional information must be 1, 2 or 3.
        % When the horizon is 1, we have a density conditional forecast. When
        % the horizon is 2, we have a soft condition with unknown central
        % tendency. When the horizon is 3, we have a soft condition with known
        % central tendency. The elements should be sorted as lower_bound,
        % central tendency, upper bound.
        % data for conditional forecast and real-time estimation
        'cond_data_ct',rise_time_series
        'cond_data_lb',rise_time_series
        'cond_data_ub',rise_time_series
        % this should be a property of solve... It should be replaced with
        % solve_exp_order. order will be reserved for higher order
        % approximation in which case it will be called solve_order
        'order',expansion_order 
        'shock_properties',shock_properties
		%
        % ====== solve options ======
        'qz_criterium',1.000001
        'dsge_varlag',4
        'dsge_var_constant',true
        % should be a direct default if the functions where it is used
        'Penalty',1e+8 
        % ====== optimization options ======
        'derivatives','symbolic' %['symbolic','numerical','automatic']
        % ====== Filtering, priors, likelihood options ======
        'hessian','fd' % (finite differences) alternatives: 'opg' (outer-product gradient)
        % ====== graphic and debugging options ======
        'graphics',[rows,cols]
        'verbose',false
        'discretize',20 % For the computation of check plots and priors
        'debug',false
        };
    MainOptions = cell2struct(MainOptions(:,2),MainOptions(:,1),1);
    obj.options=mergestructures(myoptions,MainOptions);
    
    % now update things with varargin if necessary
    obj=set_options(obj,varargin{:});
else
    nargs=length(varargin);
    if rem(nargs,2)~=0
        error([mfilename,':: arguments should enter by pairs'])
    end
    update_random_number_stream=false;
    refresh_functions=false;
    if nargs>0
        All_properties=fieldnames(obj.options);
        for ii=1:nargs/2
            propname=varargin{2*ii-1};
            if isnan(locate_variables(propname,All_properties,true))
                error([mfilename,...
                    ':: ',propname,' is not a valid option of class ',class(obj)])
            end
            propval=varargin{2*ii};
            
            if strcmp(propname,'rise_functions2disk')
                refresh_functions=true;
            end
            
            if strcmp(propname,'optimset')
                subfields=fieldnames(propval);
                for sbf=1:numel(subfields)
                    field=subfields{sbf};
                    if isfield(obj.options.(propname),field)
                        obj.options.(propname).(field)=propval.(field);
                    else
                        error([mfilename,':: ',field,' not an appropriate property of optimset '])
                    end
                end
                continue
            end
            
            if strcmp(propname,'simul_algo')
                if ~ismember(propval,...
                        {'mt19937ar','mcg16807','mlfg6331_64','mrg32k3a','shr3cong','swb2712'})
                    error([mfilename,':: simul_algo must be a member of ''mt19937ar'',''mcg16807'',''mlfg6331_64'',''mrg32k3a'',''shr3cong'',''swb2712'''])
                end
                update_random_number_stream=true;
            end
            
            if strcmp(propname,'simul_seed')
                if ~isnumeric(propval)
                    error([mfilename,':: simul_seed must must be a numerical scalar'])
                end
                update_random_number_stream=true;
            end
            
            if strcmp(propname,'estim_general_restrictions')&& ~isa(propval,'function_handle')
                if isempty(propval)
                    continue
                end
                propval=str2func(propval);
            end
            
            if strcmp(propname,'data')|| ...
                    strcmp(propname,'cond_data_ct')||...
                    strcmp(propname,'cond_data_lb')||...
                    strcmp(propname,'cond_data_ub')
                % could also check that the observable variables are included
                % in data... I do that during estimation, I guess but I can do
                % that upfront too. Don't know what is best.
                assert(isa(propval,'rise_time_series'),'data input must be a rise_time_series object')
            end
            
            if strcmp(propname,'order')
                obj.options.shock_properties=struct('name',{obj.varexo.name}',...
                    'StandardDeviation',num2cell(nan(obj.NumberOfExogenous,1)),...
                    'horizon',num2cell(propval*ones(obj.NumberOfExogenous,1)));
            end
            if strcmp(propname,'shock_properties')
                assert(numel(fieldnames(propval))==3,'shock_properties should have fields name, horizon and StandardDeviation only')
                assert(isfield(propval,'name'),'shock_properties should have a field called name')
                assert(isfield(propval,'horizon'),'shock_properties should have a field called horizon')
                nshocks=max(size(propval));
                SHOCKS={obj.options.shock_properties.name};
                for jj=1:nshocks
                    % for each one, assert that the shock exists
                    shock=propval(jj).name;
                    loc=find(strcmp(shock,{obj.varexo.name}),1);
                    if isempty(loc)
                        error([mfilename,':: ',shock,' is not a declared shock in the model file'])
                    end
                    % assign the properties
                    loc=strcmp(shock,SHOCKS);
                    obj.options.shock_properties(loc).horizon=propval(jj).horizon;
                    % update the order if necessary
                    obj.options.order=max(propval(jj).horizon,obj.options.order);
                end
                continue
            end
            
            % if applicable, assign the properties
            obj.options.(propname)=propval;
        end
    end
    
end

if update_random_number_stream
    % random number generation control
    s = RandStream.create(obj.options.simul_algo,'seed',obj.options.simul_seed);
    stream_methods=methods('RandStream');
    if ismember('setGlobalStream',stream_methods)
        RandStream.setGlobalStream(s);
    else
        RandStream.setDefaultStream(s); %#ok<SETRS>
    end
end

if refresh_functions
    obj=load_functions(obj);
end
% if ~isempty(obj.options.harvey_scale_factor)
%     obj.options.ergodic=false;
% end
