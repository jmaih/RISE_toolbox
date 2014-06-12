function [pp_i,def_i,retcode]=ergodic_parameters(Q,def,pp)
ergodic=true;
[PAI00,retcode]=initial_markov_distribution(Q,ergodic);
pp_i=0;
def_i=[];
if isempty(def{1})
    def_i=0;
end
for ii=1:size(pp,2)
    pp_i=pp_i+PAI00(ii)*pp(:,ii);
    if ~isempty(def{ii})
        def_i=def_i+PAI00(ii)*def{ii};
    end
end
end