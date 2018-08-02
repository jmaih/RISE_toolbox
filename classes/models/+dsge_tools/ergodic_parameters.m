function [pp_i,def_i,retcode]=ergodic_parameters(Q,def,pp)
% INTERNAL FUNCTION
%

ergodic=true;

[PAI00,retcode]=initial_markov_distribution(Q,ergodic);

pp_i=0;

def_i=[];

do_defs= ~isempty(def{1});

if do_defs

    def_i=0;

end

for ii=1:size(pp,2)

    pp_i=pp_i+PAI00(ii)*pp(:,ii);

    if do_defs

        def_i=def_i+PAI00(ii)*def{ii};

    end

end

end