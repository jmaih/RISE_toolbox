function nmc=initialize_markov_chain(name,nstates,endo_flag)
if nargin<3
    endo_flag=[];
end
if isempty(endo_flag)
    endo_flag=nan;
elseif ~islogical(endo_flag)
    error('endo_flag must be logical')
end
state_names=parser.create_state_list(name,nstates);
nmc=struct('name',name,...
    'number_of_states',nstates,...
    'is_endogenous',endo_flag,...
    'state_names',{state_names(:).'},...
    'state_tex_names',{state_names(:).'});
end
