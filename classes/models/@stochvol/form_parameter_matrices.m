function [A,RHO_A,THETA_A,...
    SIG,RHO_SIG,THETA_SIG,...
    OMG,RHO_OMG,THETA_OMG]=form_parameter_matrices(obj)
obj=solve(obj);

SIG=obj.solution.sig;
RHO_SIG={};
THETA_SIG={};
if obj.time_varying_parameters(2)
    RHO_SIG=obj.solution.rho_sig;
    THETA_SIG=obj.solution.theta_sig;
end

OMG=obj.solution.omg;
RHO_OMG={};
THETA_OMG={};
if obj.time_varying_parameters(3)
    RHO_OMG=obj.solution.rho_omg;
    THETA_OMG=obj.solution.theta_omg;
end

A=cell(1,obj.markov_chains.regimes_number);
RHO_A=A;
THETA_A=A;
for ireg=1:obj.markov_chains.regimes_number
    for ilag=1:obj.nlags
        A{ireg}=[A{ireg},obj.solution.(sprintf('a%0.0f',ilag)){ireg}];
        if obj.time_varying_parameters(1)
            RHO_A{ireg}=[RHO_A{ireg},obj.solution.(sprintf('rho_a%0.0f',ilag)){ireg}];
            THETA_A{ireg}=[THETA_A{ireg},obj.solution.(sprintf('theta_a%0.0f',ilag)){ireg}];
        end
    end
    if isfield(obj.solution,'c')
        % deterministic terms
        A{ireg}=[A{ireg},obj.solution.c{ireg}];
        if obj.time_varying_parameters(1)
            RHO_A{ireg}=[RHO_A{ireg},obj.solution.rho_c{ireg}];
            THETA_A{ireg}=[THETA_A{ireg},obj.solution.theta_c{ireg}];
        end
    end
end

end