allItems=struct('name',{},...
    'type',{},...
    'observed',{},...
    'switching',{},...
    'tex',{},...
    'log_var',{},...
    'switch',{},...
    'id',{},...
    'markov_chain',{});

iter=0;
jter=0;
obs_iter=0;
    type=1;
for ivar=1:numel(endo_list)
    iter=iter+1;
    jter=jter+1;
    allItems(iter).name=pname{ivar};
    allItems(iter).tex=ptexname{ivar};
    allItems(iter).type=type;
    allItems(iter).switch=nan;
    allItems(iter).id=jter;
    if isobserved()
        obs_iter=obs_iter+1;
        allItems(iter).id=allItems(iter).id+obs_iter*1i;
    end
    allItems(iter).markov_chain=0;
end

for ivar=1:numel(exo_list)
end

for ivar=1:numel(param_list)
end