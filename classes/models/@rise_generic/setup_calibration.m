function obj=setup_calibration(obj,Calibration)
% Calibration is either a struct or a cell of the form {names,paramvector}
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
    pnames=Calibration{1};
    % make sure the names are not in the tex form
    pnames=parser.param_name_to_valid_param_name(pnames);
    param_draw=Calibration{2};
    Calibration=struct();
    for iname=1:numel(pnames)
        Calibration.(pnames{iname})=param_draw(iname);
    end
end

pnames=fieldnames(Calibration);

% push the calibration
%---------------------
[position,regime_states]=generic_tools.parameter_position_and_regimes(pnames,...
    param_names,governing_chain,chain_names,grand_chains_to_small,regimes);

fields=fieldnames(Calibration);
for ii=1:numel(fields)
    tmp=Calibration.(fields{ii});
    obj.parameter_values(position(ii),regime_states{ii})=tmp;
end
end