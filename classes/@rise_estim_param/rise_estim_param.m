classdef rise_estim_param < rise_param
    properties
        plb
        pub
        distribution
        interval_probability
    end
    properties(Dependent)
    end
    properties (SetAccess = private, Hidden = true)
        ab_space
    end
    properties(SetAccess=protected)
        a
        b
        c=nan;
        d=nan;
        lb=[];
        ub=[];
        mode=[];
        mode_std=[];
        mean=[];
        quantiles=[];
        post_sim_mode=[];
        prior_mean
        prior_standard_deviation
    end
    methods
        % set property utility
        obj=set_properties(obj,varargin)
        plot(obj)
        % constructor
        function obj=rise_estim_param(name,tex_name,id,value,plb,pub,distrib,prob,prior_trunc,c,d)
			obj=obj@rise_param('name',name,'tex_name',tex_name,'id',id,'startval',value);
            obj.plb=plb;
            obj.pub=pub;
            obj.distribution=strrep(distrib,'_pdf','');
            obj.interval_probability=prob;
            if prob<=0||prob>1
                error([mfilename,':: probability for parameter ',name,' should be in (0,1]'])
            end
            if nargin>9
                obj.c=c;
                if nargin>10
                    obj.d=d;
                end
            end
            %=========================
            % for truncation
            invgamma_upper_bound_truncation=10;
            %=========================
            % find the hyperparameters
            [obj.a,... 
                obj.b,... 
                moments,... %
                ffinal,...
                obj.ab_space]=distributions.(obj.distribution)(plb,pub,prob,obj.c,obj.d); 
            obj.prior_mean=moments.mean;
            obj.prior_standard_deviation=moments.sd;
            disp([' parameter: ',upper(name),', density:',upper(obj.distribution),...
                ', hyperparameters: [',num2str(obj.a),' ',num2str(obj.b),'],',...
                'convergence ',num2str(ffinal)])
            % get the functions of the distribution
            [~,~,icdfn]=distributions.(obj.distribution)();
            bounds=[icdfn(prior_trunc,obj.a,obj.b,obj.c,obj.d),...
                icdfn(1-prior_trunc,obj.a,obj.b,obj.c,obj.d)];
            if isempty(obj.prior_mean)||isempty(obj.prior_standard_deviation)
                disp([mfilename,'(GENTLE WARNING):: for these hyperparameters, the distribution ',...
                        'does not have well-defined moments'])
            end
            if ismember(obj.distribution,{'inv_gamma'})
                if bounds(2)>invgamma_upper_bound_truncation
                    the_message=[mfilename,'(GENTLE WARNING):: upper bound of inverse gamma distribution ',...
                        'truncated at ',num2str(invgamma_upper_bound_truncation)];
                end
                bounds(2) = min(bounds(2),invgamma_upper_bound_truncation);
            else
                the_message='';
            end
            % if the distribution has been truncated, say it here.
            disp(the_message)
            % check that the starting value is not outside the bounds
            if any(obj.startval<bounds(1))||any(obj.startval>bounds(2))
                error([mfilename,':: parameter ',name,' (',num2str(obj.startval),') outside its bounds [',num2str(bounds),']'])
            end
            obj.lb = bounds(1);
            obj.ub = bounds(2);
        end
        function obj=reset_start_value(obj,name,value,id)
            % this function resets the start values for estimation. the id
            % should be provided if there are several objects with the same
            % name in the vector of 'rise_estim_param' objects. The ids
            % are the same as the order of the estimated parameters in your
            % model file. This is formalized in private function
            % rise/format_parameters
            try
                narginchk(3,4)
            catch me
                % for backward compatibility
                warning(me.message)
                error(nargchk(3,4,nargin,'string')) %#ok<NCHKN>
            end
            if nargin<4
                id=[];
            end
            par_names={obj.name};
            if ischar(name)
                name=cellstr(name);
            end
            bad_names=~ismember(name,par_names);
            if any(bad_names)
                disp(name(bad_names))
                error([mfilename,':: the parameter names listed above are not in the list of estimated parameters'])
            end
            if numel(name)~=numel(value)
                error([mfilename,':: number of names should be equal to number of values'])
            end
            ids=[obj.id];
            if ~isempty(id)
                bad_ids=~ismember(id,ids);
                if any(bad_ids)
                    disp(id(bad_ids))
                    error([mfilename,':: the ids above are not in the range of valid ids'])
                end
            end
            for ii=1:numel(name)
                name_loc=strcmp(name{ii},par_names);
                if numel(name_loc)>1
                    if isempty(id)
                        error([mfilename,':: multiple occurrences of ''',name{ii},''' in the list of estimated parameters. ids should be provided'])
                    end
                    id_i=id(ii);
                else
                    id_i=find(name_loc);
                end
                if value(ii)>obj(id_i).ub||value(ii)<obj(id_i).lb
                    error([mfilename,':: value for parameter ''',name{ii},''' outside its bounds ([',num2str(obj(id_i).lb),' ',num2str(obj(id_i).ub),'])'])
                end
                obj(id_i).startval=value(ii);
            end
        end
        function obj=reset_bounds(obj,low,high)
            if nargin<3
                high=nan;
                if nargin<2
                    low=nan;
                end
            elseif nargin>3
                error([mfilename,':: number of arguments cannot exceed 3'])
            end
            set_input('lb',low)
            set_input('ub',high)
            function set_input(where,x)
                if isnan(x)||isinf(x)
                    warning([mfilename,':: nan or inf new value not set for ',where])
                elseif ~isreal(x)
                    error([mfilename,':: input must be real'])
                else
                    nobj=numel(obj);
                    for iobj=1:nobj
                        switch where
                            case 'lb'
                                if x>obj(iobj).ub
                                    error([mfilename,':: attempting to set a value for the lower bound that exceeds the upper bound'])
                                end
                            case 'ub'
                                if x<obj(iobj).lb
                                    error([mfilename,':: attempting to set a value for the upper bound that exceeds the lower bound'])
                                end
                        end
                        obj(iobj).(where)=x;
                    end
                end
            end
        end
    end
end


