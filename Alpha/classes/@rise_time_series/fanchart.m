function out=fanchart(data,ci)
out=struct();
nvar=data.NumberOfVariables;
if nvar>1
    vnames=data.varnames;
    for ivar=1:nvar
        out.(vnames{ivar})=fanchart(window(data,[],[],vnames{ivar}),ci);
    end
    return
end
ci=sort(ci(:));
datax=squeeze(double(data));
probs = 0.5*(1+[-ci,ci]);
probs=[flipud(probs(:,1));probs(:,2)];
emp_moms=distributions.empirical_moments(datax,[],[],probs);
out=struct();
out.ci=ci;
out.mean=[emp_moms.mean];
out.median=[emp_moms.median];
out.variance=[emp_moms.variance];
out.quantiles=vertcat(emp_moms.quantiles);
out.prob_index=emp_moms(1).prob_index;
out.probs=emp_moms(1).probs;
out.x=data.date_number;