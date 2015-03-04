function obj=set(obj,varargin)
% set - sets options for RISE models
%
% Syntax
% -------
% ::
%   
%   obj=set(obj,varargin)
%
% Inputs
% -------
%
% - **obj** [rise|dsge|svar|rfvar]: model object
%
% - **varargin** : valid input arguments coming in pairs.
%
% Outputs
% --------
%
% - **obj** [rise|dsge|svar|rfvar]: model object
%
% More About
% ------------
%
% - one can force a new field into the options by prefixing it with a '+'
%   sign. Let's say yourfield is not part of the options and you would like
%   to force it to be in the options because it is going to be used in some
%   function or algorithm down the road. Then you can run
%   m=set(m,'+yourfield',value). then m will be part of the new options.
%
% Examples
% ---------
%
% See also: 

if isempty(obj)
    obj=struct();
    return
end
nobj=numel(obj);
if nobj>1
    is_legend=false;
    leg_pos=[];
    for ii=1:2:length(varargin)
        new_legend=strcmp(varargin{ii},'legend');
        if new_legend
            if is_legend
                error('more than one legends found in same call')
            end
            is_legend=new_legend;
            leg_pos=ii+1;
            if ~iscellstr(varargin{leg_pos})
                error('with multiple objects, legend must be a cell array of strings')
            end
            if numel(unique(varargin{leg_pos}))~=numel(varargin{leg_pos})
                error('legends must be unique to each model')
            end
        end
    end
    cell_argin=varargin;
    for iobj=1:nobj
        if is_legend
            iobj_leg=varargin{leg_pos}{iobj};
            if isempty(iobj_leg)
                iobj_leg=sprintf('model_%0.0f',iobj);
            end
            cell_argin{leg_pos}=iobj_leg;
        end
        obj(iobj)=set(obj(iobj),cell_argin{:});
    end
    return
end

initialize=length(varargin)==1 && strcmp(varargin{1},'initialize');
if initialize
    varargin=[];
end

if ~isempty(varargin) && isstruct(varargin{1})
    % convert to cell
    tmp=struct2cell(varargin{1});
    fnames=fieldnames(varargin{1});
    varargin=transpose([fnames(:),tmp(:)]);
    varargin=varargin(:)';
end

nn=length(varargin);
if rem(nn,2)
    error('arguments must come in pairs')
end

update_random_number_stream=false;
if isempty(obj.options)
    set_all_options()
end

if nn==0 && ~initialize
    % list all the settable options
    %------------------------------
elseif nn==1
    % list the options for the particular field
    %------------------------------------------
    error('something should be implemented here. Please send a reminder to junior.maih@gmail.com')
else
    options_fields=fieldnames(obj.options);
    for ii=1:2:nn
        propname=varargin{ii};
        propval=varargin{ii+1};
        forced=strcmp(propname(1),'+');
        if forced
            propname=propname(2:end);
        end
        if ~isvarname(propname)
            error([propname,' is not a valid variable name'])
        end
        set_one_at_a_time();
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
    function set_one_at_a_time()
        if strcmp(propname,'filters')
            obj.filtering=propval;
        elseif strcmp(propname,'tex_name')
            set_tex_names()
        elseif any(strcmpi(propname,{'priors','estim_priors'}))
            obj=setup_priors(obj,propval);
        elseif strcmpi(propname,'parameters')
            % propval is either a struct or a cell of the form {names,paramvector}
            obj=setup_calibration(obj,propval);
        elseif strcmpi(propname,'legend')
            if ~ischar(propval)
                error('legend must be a string')
            end
            obj.legend=propval;
        elseif any(strcmp(propname,options_fields))||forced
            set_one_option()
        else
            error(['',propname,''' is not a settable option or property of class ',class(obj)])
        end
        function set_tex_names()
            if ~isstruct(propval)
                if ~iscellstr(propval) && size(propval,2)==2
                    error('propval must be a structure of a cellstr with two columns')
                end
                propval=cell2struct(propval(:,2),propval(:,1),1);
            end
            fields=fieldnames(propval);
            for ifield=1:numel(fields)
                thisfield=fields{ifield};
                loc=find(strcmp(thisfield,obj.endogenous.name));
                if isempty(loc)
                    loc=find(strcmp(thisfield,obj.exogenous.name));
                    if isempty(loc)
                        loc=find(strcmp(thisfield,obj.parameters.name));
                        if isempty(loc)
                            loc=find(strcmp(thisfield,obj.markov_chains.chain_names));
                            if isempty(loc)
                                loc=find(strcmp(thisfield,obj.markov_chains.state_names));
                                if isempty(loc)
                                    error(['"',thisfield,'" is not recognized as an endog, an exog, a param, a chain name or a regime name'])
                                end
                                obj.markov_chains.state_tex_names{loc}=propval.(thisfield);
                                continue
                            else
                                obj.markov_chains.chain_tex_names{loc}=propval.(thisfield);
                                continue
                            end
                        else
                            type='parameters';
                        end
                    else
                        type='exogenous';
                    end
                else
                    type='endogenous';
                end
                obj.(type).tex_name{loc}=propval.(thisfield);
                if ~isempty(obj.observables.name) && ...
                        any(strcmp(type,{'endogenous','exogenous'}))
                    loc=find(strcmp(thisfield,obj.observables.name));
                    if ~isempty(loc)
                        obj.observables.tex_name{loc}=propval.(thisfield);
                    end
                end
            end
        end
        function set_one_option()
            absorbed=false;
            if strcmp(propname,'optimset')
                subfields=fieldnames(propval);
                for sbf=1:numel(subfields)
                    field=subfields{sbf};
                    if isfield(obj.options.(propname),field)
                        obj.options.(propname).(field)=propval.(field);
                        absorbed=true;
                    else
                        error([mfilename,':: ',field,' not an appropriate propname of optimset '])
                    end
                end
            elseif any(strcmp(propname,{'estim_start_date','estim_end_date',...
                    'data','data_demean','kf_presample'}))
                obj.data_are_loaded=false;
            elseif any(strcmp(propname,{'data','data_cond_ct','data_cond_lb','data_cond_ub'}))
                % could also check that the observable variables are included
                % in data... I do that during estimation, I guess but I can do
                % that upfront too. Don't know what is best.
                propval=ts.collect(propval);
            elseif strcmp(propname,'simul_algo')
                if ~ismember(propval,...
                        {'mt19937ar','mcg16807','mlfg6331_64','mrg32k3a','shr3cong','swb2712'})
                    error([mfilename,':: simul_algo must be a member of ''mt19937ar'',''mcg16807'',''mlfg6331_64'',''mrg32k3a'',''shr3cong'',''swb2712'''])
                end
                update_random_number_stream=true;
            elseif strcmp(propname,'simul_seed')
                if ~isnumeric(propval)
                    error([mfilename,':: simul_seed must must be a numerical scalar'])
                end
                update_random_number_stream=true;
            elseif strcmp(propname,'estim_general_restrictions')&& ~isa(propval,'function_handle')
                if isempty(propval)
                    absorbed=true;
                else
                    if ischar(propval)
                        propval=str2func(propval);
                    elseif iscell(propval) && ischar(propval{1})
                        propval{1}=str2func(propval{1});
                    end
                end
            end

            if ~absorbed
                obj.options.(propname)=propval;
            end
        end
    end
    function set_all_options()
        % create a dummy object to read the options from various
        % objects, functions and sub-functions
        obj_class=class(obj);
        dum_obj=eval([obj_class,'.empty(0,1)']);
        method_list=methods(dum_obj);
        restricted_list={'Contents',class(obj),'horzcat','cat','vertcat',...
            'template'};%,'prior_plots','initialize_solution_or_structure'
        method_list(ismember(method_list,restricted_list))=[];
        % add some private methods to the list
        method_list=[method_list(:)','load_data'];%;'load_functions'
        myoptions=struct();
        for imeth=1:numel(method_list)
            myoptions=utils.miscellaneous.mergestructures(myoptions,...
                dum_obj.(method_list{imeth})...
                );
        end
        
        % Create the main options
        update_random_number_stream=true;
        % shocks are allowed to be anticipated.
        rows=4;cols=3;
        if isa(obj,'dsge') % isfield(obj,'filename') does not work!!!
            results_folder=obj.filename;
        else
            ymdhms=datevec(now);
            year=sprintf('%0.0f',ymdhms(1));
            
            month=sprintf('%0.0f',ymdhms(2));
            if length(month)==1,month=['0',month];end
            
            day=sprintf('%0.0f',ymdhms(3));
            if length(day)==1,day=['0',day];end
            
            hour=sprintf('%0.0f',ymdhms(4));
            if length(hour)==1,hour=['0',hour];end
            
            minutes=sprintf('%0.0f',ymdhms(5));
            if length(minutes)==1,minutes=['0',minutes];end
            results_folder=[year,month,day,hour,minutes];
        end
        MainOptions={
            % set output folder name
            'results_folder',results_folder
            % ====== solve options ======
            'dsge_varlag',4
            'dsge_var_constant',true
            % ====== graphic and debugging options ======
            'graphics',[rows,cols]
            'verbose',false
            'debug',false
            };
        MainOptions = cell2struct(MainOptions(:,2),MainOptions(:,1),1);
        obj.options=utils.miscellaneous.mergestructures(myoptions,MainOptions);
    end
end

