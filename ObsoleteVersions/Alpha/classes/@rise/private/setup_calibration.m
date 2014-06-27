function obj=setup_calibration(obj,Calibration)
% Calibration is either a struct or a cell of the form {names,paramvector}
param_names=obj.parameters.name;

par_nbr=sum(obj.parameters.number);
regimes=cell2mat(obj.markov_chains.regimes(2:end,2:end));
reg_nbr=size(regimes,1);
if isempty(obj.parameter_values)
    obj.parameter_values=nan(par_nbr,reg_nbr);
end
chain_names=obj.markov_chains.chain_names;

% Transform to struct if necessary
%---------------------------------
if iscell(Calibration)
    pnames=Calibration{1};
    % make sure the names are not in the tex form
    pnames=parser.param_name_to_valid_param_name(pnames);
    param_draw=Calibration{2};
    Calibration=struct();
    for iname=1:numel(pnames)
        Calibration.(pnames{iname})=param_draw(iname);
    end
end

% push the calibration
%---------------------
fields=fieldnames(Calibration);
for ii=1:numel(fields)
    tmp=Calibration.(fields{ii});
    [position,regime_states]=position_and_regimes(fields{ii});
    obj.parameter_values(position,regime_states)=tmp;
end

    function [position,regime_states]=position_and_regimes(pname)
        position=find(strcmp(pname,param_names));
        if ~isempty(position)
            state=1;
            chain='const';
            chain_id=find(strcmp(chain,chain_names));
        else
            ptex=parser.valid_param_name_to_tex_name(pname,chain_names);
            left_par=strfind(ptex,'(');
            right_par=strfind(ptex,')');
            pname=ptex(1:left_par-1);
            comma=strfind(ptex,',');
            position=find(strcmp(pname,param_names));
            if isempty(left_par)||isempty(right_par)||isempty(comma)||isempty(position)
                error(['"',ptex(1:left_par-1),'" is not recognized as a parameter name'])
            end
            chain=ptex(left_par+1:comma-1);
            state=str2double(ptex(comma+1:right_par-1));
            chain_id=find(strcmp(chain,chain_names));
        end
        govChain=obj.parameters.governing_chain(position);
        if ~(chain_id==govChain)
            error(['parameter ',pname,' is not controlled by ',chain_names{chain_id}])
        end
        % locate the state in the regimes
        regime_states=find(regimes(:,chain_id)==state);
        if isempty(regime_states)
            error([sprintf('%0.0f',state),' is not a valid state for parameter ',pname])
        end
    end
end