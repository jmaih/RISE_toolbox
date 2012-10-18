function [p,fval,HESS]=studga_gate_iris(Objective,p0,pl,pu,optim_set)

fields={'colony_size','max_iter','max_time','max_fcount'};
fields_equiv={'MaxNodes','MaxIter','MaxTime','MaxFunEvals'};
fields_defaults={20,10000,60*60,inf};

for ii=1:numel(fields_equiv)
    if isfield(optim_set,fields_equiv{ii}) && ...
            ~isempty(optim_set.(fields_equiv{ii}))
        options.(fields{ii})=optim_set.(fields_equiv{ii});
    else
        options.(fields{ii})=fields_defaults{ii};
    end
end

[p,fval]=studga_gate(Objective,p0(:),[],pl(:),pu(:),options);

HESS=finite_difference_hessian(Objective,p);

p=reshape(p,size(p0));