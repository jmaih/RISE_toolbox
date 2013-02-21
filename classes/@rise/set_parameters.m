function obj=set_parameters(obj,param_struct)
% There are 5 ways to assign parameters, depending on what exactly you
% would like to do:
% Case 1: a structure with names of the parameters as fields. In case of a 
% regime switching model, entering a parameter name as given in the model 
% implies that the value assigned will be the same in all the regimes. In
% other words, this is the same as assuming that the parameter is controled
% by the "const" (constant) chain. e.g. param_struct.alpha=1; In case the
% parameter is controled by a different markov process and hence assumes
% different values across the states of the commanding chain, the name
% entered in the structure is augmented with the name of the chain and the
% state. e.g. param_struct.alpha_coef_3=1; This says that "alpha" is
% controled by a markov chain with name "coef" and takes on a value of "1"
% in state "3"
%
% Case 2: a 2-element cell where the first cell is the ids in the matrix of
%  parameters, which is of dimensions n_params x n_regimes. The second
%  element is the vector of values to be assigned. This case is also not
%  advisable and perhaps I should remove this possibility.
%
% Case 3: the model is set up for estimation. Then param_struct is a any of
% the following strings: 'startval', 'mode', 'mean' or 'post_sim_mode',
% depending on whether the model has been estimated or not.
%
% Case 4: the model is set up for estimation. Then param_struct can be a
% vector of numbers matching the estimation restrictions. This is the
% typical option used by the algorithm during estimation and is not
% advisable to use it otherwise.
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
    true_names={obj.parameters.name};
    param_struct=param_struct(:);
    npar=numel(param_struct);   
    direct_assignment()
    is_done=true;
elseif iscell(param_struct) % Case 4
    if numel(param_struct)~=2
        help('set_parameters')
        error([mfilename,':: follow one of the options above'])
    end
    id=param_struct{1};
    values=param_struct{2};
    if numel(id)~=numel(values)
        error([mfilename,':: # locations different from # values'])
    end
    if max(id)>param_nbr*obj.NumberOfRegimes
        error([mfilename,':: elements in the first cell cannot be greater than ',...
            int2str(param_nbr*obj.NumberOfRegimes)])
    end
else
    help('set_parameters')
    error([mfilename,':: follow one of the options above'])
end

if ~is_done
    % now set the parameters
		loc= strcmp('startval',obj.parameters_image(1,:));
		obj.parameters_image{2,loc}(id)=values(:);
	% now push back the parameters into the object if possible
	obj=rehash(obj);
end

    function direct_assignment()
        pnames=fieldnames(param_struct);
        default_cols=1:obj.NumberOfRegimes;
        % I should have a more permanent way to locate the chains and their
        % states that does not require going into endless depths of the
        % structures. A cell with names on top should suffice...
        chains_names={obj.markov_chains.name};
        loc= strcmp('startval',obj.parameters_image(1,:));
        for ipar=1:npar
            cols=default_cols;
            myparam=pnames{ipar};
            name_loc=find(strcmp(myparam,true_names));
            if isempty(name_loc)
                underscores=find(myparam=='_');
                if numel(underscores)<2
                    error(['no parameter with name ',myparam,' found'])
                end
                myname=myparam(1:underscores(end-1)-1);
                name_loc=find(strcmp(myname,true_names));
                if isempty(name_loc)
                    error(['no parameter with names ',myparam,' or ',myname,' found'])
                end
                chain=myparam(underscores(end-1)+1:underscores(end)-1);
                chain_loc=find(strcmp(chain,chains_names));
                if isempty(chain_loc)
                    error(['',chain,''' is not recognized as a chain name'])
                end
                state=str2double(myparam(underscores(end)+1:end));
                cols=find(obj.Regimes(:,chain_loc)==state);
                if isempty(cols)
                    error(['state ',num2str(state),' outside the state range of markov chain ',chain])
                end
            end
            obj.parameters_image{2,loc}(name_loc,cols)=param_struct.(myparam);
        end
        obj=rehash(obj);
    end
end