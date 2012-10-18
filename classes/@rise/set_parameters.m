function obj=set_parameters(obj,param_struct)
% There are 5 ways to assign parameters, depending on what exactly you
% would like to do:
% Case 1: the model is set up for estimation. Then param_struct is a any of
% the following strings: 'startval', 'mode', 'mean' or 'post_sim_mode',
% depending on whether the model has been estimated or not.
%
% Case 2: the model is set up for estimation. Then param_struct can be a
% vector of numbers matching the estimation restrictions. This is the
% typical option used by the algorithm during estimation and is not
% advisable to use it otherwise.
%
% Case 3: a structure with fields (name,value,regime) 
% example: param_struct=struct('name',{'alpha','beta'},...
%                           'regime',{1,2},...
%                           'value',{0.1,5});
%
% Case 4: a 2-element cell where the first cell is the ids in the matrix of
%  parameters, which is of dimensions n_params x n_regimes. The second
%  element is the vector of values to be assigned. This case is also not
%  advisable and perhaps I should remove this possibility.
%
% Case 5: a 3-column cell with headers 'name', 'value','regime' (order does not
% matter).
%    example: params={'name', 'value','regime'
%                     'alpha',  0.1  ,  1
%                     'beta',   5    ,  1}
%
% See also : assign_estimates

if isempty(obj)
    obj=struct();
    return
end

param_nbr=obj.NumberOfParameters;
is_done=false;
if ischar(param_struct) % Case 1
    % then param_struct is one of the following : startval, mode, mean
    params=vertcat(obj.estimated_parameters.(param_struct));
    obj=assign_estimates(obj,params);
    is_done=true;
elseif isvector(param_struct) && isnumeric(param_struct) % Case 2
    % this is used especially during estimation as restrictions are needed
    % to know the names to assign the parameters to.
    obj=assign_estimates(obj,param_struct(:));
    is_done=true;
elseif isstruct(param_struct) % Case 3
    % 1- vectorize
    param_struct=param_struct(:);
    if isfield(param_struct,'name')
        names={param_struct.name};
    else
        error([mfilename,':: argument must be a structure with ''name'' as one of the fields'])
    end
    if isfield(param_struct,'value')
        values=[param_struct.value];
    else
        error([mfilename,':: argument must be a structure with ''value'' as one of the fields'])
    end
    % number of parameters to set
    npar=numel(param_struct);
    if isfield(param_struct,'regime')
        regimes=[param_struct.regime];
    else
        regimes=ones(npar,1);
    end
    % check consistency of structure
    unique_names=unique(names);
    parlocs=locate_variables(names,{obj.parameters.name});
    for ii=1:numel(unique_names)
        name_id=find(strcmp(unique_names{ii},names));
        if numel(name_id)>1
            reg_i=regimes(name_id);
            if numel(reg_i)~=numel(unique(reg_i))
                error([mfilename,':: parameter ',unique_names{ii},' is listed at least twice in the same regime'])
            end
            % then the parameter is switching. We have to confirm it here
            % in case the model was originally not parameterized... It is
            % redundant if the model was originally parameterized and slows
            % things down a bit under estimation.
            if ~obj.estimation_under_way
                id=parlocs(name_id(1)); % <--- id=parlocs(ii);
                name_=unique_names{ii}; % same as obj.parameters(id).name;
                tex_name_=obj.parameters(id).tex_name;
                startval_=obj.parameters(id).startval;
                obj.parameters(id)=rise_param('name',name_,'tex_name',tex_name_,'id',id,'startval',startval_,'is_switching',true);
            end
        end
    end
    if max(regimes)>obj.NumberOfRegimes
        error([mfilename,':: the largest regime in your parameters (',...
            int2str(max(regimes)),') exceeds the number of regimes ',int2str(obj.NumberOfRegimes)])
    end
    id=parlocs(:)+(regimes(:)-1)*param_nbr;
elseif iscell(param_struct) % Case 4
    if numel(param_struct)==2
        id=param_struct{1};
        values=param_struct{2};
        if numel(id)~=numel(values)
            error([mfilename,':: # locations different from # values'])
        end
        if max(id)>param_nbr*obj.NumberOfRegimes
            error([mfilename,':: elements in the first cell cannot be greater than ',...
                int2str(param_nbr*obj.NumberOfRegimes)])
        end
    elseif size(param_struct,2)==3 && all(ismember({'name','value','regime'},param_struct(1,:))) % Case 5
        param_struct=cell2struct(param_struct(2:end,:),param_struct(1,:),2);
        obj=set_parameters(obj,param_struct);
        is_done=true;
    else
        error([mfilename,':: wrong specification of the parameter input'])
    end
else
    help('set_parameters')
    error([mfilename,':: follow one of the options above'])
end

if ~is_done
    % now set the parameters
    if obj.estimation_under_way
		loc=find(strcmp('startval',obj.parameters(1,:)));
		obj.parameters{2,loc}(id)=values(:);
	else
	    Mparams=vertcat(obj.parameters.startval);
	    Mparams(id)=values(:);
        % pushing the parameters in this way is very expensive... I need to
        % find a better way
        Mparams=mat2cell(Mparams,ones(param_nbr,1),obj.NumberOfRegimes);
        [obj.parameters(:).startval]=deal(Mparams{:});
    end
end

end

