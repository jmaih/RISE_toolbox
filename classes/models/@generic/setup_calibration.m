function obj=setup_calibration(obj,Calibration)
% setup_calibration -- set parameters
%
% ::
%
%
%   obj=setup_calibration(obj,Calibration)
%
% Args:
%
%    - **obj** [generic]: model object
%
%    - **Calibration** [struct|cell]: calibration to push. There are two
%    possibilities:
%      - Calibration is a struct: the fields are the parameter names and each
%      name contains a parameter value.
%      - Calibration is a cell: the cell has two columns. The first column
%      holds the names of the parameters to change and the second column their
%      values.
%
% Returns:
%    :
%
%    - **obj** [generic]: model object
%
% Note:
%
%    - the parameter values could be vectors if e.g. we want to take the entire
%    parameterization of one model and push it into an identical model. But it
%    is not allowed to have situations where one parameter is a vector and
%    some other is not. Then RISE will complain that the parameter is not
%    controlled by the const markov chain.
%
% Example:
%
%    See also:

if isempty(Calibration)
    return
end
param_names=obj.parameters.name;

grand_chains_to_small=obj.markov_chains.grand_chains_to_small;
par_nbr=sum(obj.parameters.number);
regimes=cell2mat(obj.markov_chains.small_markov_chain_info.regimes(2:end,2:end));
reg_nbr=size(regimes,1);
if isempty(obj.parameter_values)
    obj.parameter_values=nan(par_nbr,reg_nbr);
end
chain_names=obj.markov_chains.small_markov_chain_info.chain_names;
governing_chain=obj.parameters.governing_chain;
% Transform to struct if necessary
%---------------------------------
if iscell(Calibration)
    pnames=Calibration(:,1);
    if iscell(pnames{1})
        pnames=pnames{1};
    end
    % make sure the names are not in the tex form
    pnames=parser.param_texname_to_param_name(pnames);
    param_draw=Calibration(:,2);
    if numel(param_draw)==1 && numel(param_draw{1})>1
        param_draw=num2cell(param_draw{1});
    end
    Calibration=struct();
    for iname=1:numel(pnames)
        Calibration.(pnames{iname})=param_draw{iname,:};
    end
end

pnames=fieldnames(Calibration);

if vector_style()
    % do nothing, the parameters are pushed automagically
else
    % find positions
    %---------------
    [position,regime_states]=generic_tools.parameter_position_and_regimes(pnames,...
        param_names,governing_chain,chain_names,grand_chains_to_small,regimes);
    
    % push the calibration
    %---------------------
    for ii=1:numel(pnames)
        tmp=Calibration.(pnames{ii});
        obj.parameter_values(position(ii),regime_states{ii})=tmp;
    end
end

    function flag=vector_style()
        silent=true;
        position=locate_variables(pnames,param_names,silent);
        flag=~any(isnan(position));
        if flag
            regimes_number=obj.markov_chains.regimes_number;
            parameter_values=obj.parameter_values;
            for ip=1:numel(pnames)
                values=Calibration.(pnames{ip});
                if ~(size(values,2)==regimes_number);
                    parameter_values=[];
                    break
                end
                parameter_values(position(ip),:)=values;
            end
            flag=~isempty(parameter_values);
            if flag
                obj.parameter_values=parameter_values;
            end
        end
    end
end