function [mysimul,retcode]=resimulate(m,f,x0h,fxlogvars)

if nargin<4
    
    fxlogvars=[];
    
end

if isempty(fxlogvars)

    fxlogvars = m.endogenous.is_log_var;

end

if nargin<3||isempty(x0h)
    
    x0h=0*permute(dsge_tools.get.stst_(m),[1,3,2]);
    
end

npar=dsge_tools.get.number_of_parameterizations(m);

mysimul=cell(1,npar);

retcode=zeros(1,npar);

for ip=1:npar

    [mysimul{ip},retcode(ip)]=engine(dsge_tools.solve.load_state_space(m,ip),...
        f,x0h,fxlogvars);

end

mysimul=ts.format_multiple_outputs(mysimul);

end

function [mysimul,retcode]=engine(m,f,x0h,fxlogvars)

tol=1e-4;

mysimul=struct();

retcode = 0;

[T,~,stst,state_vars_location]=dsge_tools.solve.load_growth_solution(m);

ff=dsge_tools.filter.histdec.forecasting_function(T,stst,state_vars_location);

endo_list=m.endogenous.name;

endo_nbr=numel(endo_list);

exo_list=m.exogenous.name;

endo_list0=fieldnames(f.smoothed_variables);
% some variables in the old model may be absent from the present model
ping=ismember(endo_list0,endo_list);

h=m.markov_chains.regimes_number;

nsmpl=get(f.smoothed_shocks.(exo_list{1}),'NumberOfObservations');

[DataEndo]=dsge_tools.filter.collect_regime_specific_data(f.smoothed_variables,endo_list0,h,nsmpl);

endo_list0=endo_list0(ping);

if ~isempty(fxlogvars)
    
    DataEndo(fxlogvars,:,:)=log(DataEndo(fxlogvars,:,:));
    
end
    
oldEndoLocs=locate_variables(endo_list0,endo_list);

x0h(oldEndoLocs,1,:)=DataEndo(ping,1,:);

delta=zeros(endo_nbr,nsmpl,h);

if isfield(f,'smoothed_switch_impact')
    
    delta(oldEndoLocs,:,:)=dsge_tools.filter.collect_regime_specific_data(f.smoothed_switch_impact,endo_list0,h,nsmpl);
    
end

shocks=dsge_tools.filter.collect_regime_specific_data(f.smoothed_shocks,exo_list,h,nsmpl);

Probs=dsge_tools.filter.collect_regime_specific_data(f.smoothed_regime_probabilities,...
    fieldnames(f.smoothed_regime_probabilities),1,nsmpl);

constraints=[];

if ~isempty(m.options.simul_constraints)
    
    simul_constraints=m.options.simul_constraints.parsed; 
    
    [simul_constraints]=obcas.memoize_constants(simul_constraints,...
        m.constants.values);
    
    h=size(Probs,1);
    
    iscontrained_regime=false(1,h);
    
    for iv=simul_constraints.constr_vars
        
        for ih=1:h
            
            if iscontrained_regime(ih),continue,end
            
            d=stst(iv,ih)-simul_constraints.recovfun(ih);
            
            iscontrained_regime(ih)=abs(d)<tol;
            
        end
        
    end
    
    constraints={simul_constraints,iscontrained_regime};
    
end


[x,Probs]=dsge_tools.filter.histdec.simulate(ff,x0h,shocks,Probs,delta,...
    constraints);

xx=dsge_tools.filter.expected_value(x,Probs,endo_nbr,h);

if any(isnan(xx(:)))||~all(isfinite(xx(:)))
    
    retcode=7004; 
    
    return
    
end

start_date=f.smoothed_regime_probabilities.regime_1.firstDateObj;

tmplt=ts(start_date,zeros(nsmpl,1));

xx(fxlogvars,:)=exp(xx(fxlogvars,:));

for iii=1:endo_nbr
    
    v=endo_list{iii};
    
    mysimul.(v)=set(tmplt,'data',xx(iii,:).');
    
end

end