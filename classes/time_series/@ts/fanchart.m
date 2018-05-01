function out=fanchart(this,ci)
% H1 line
%
% ::
%
%
% Args:
%
% Returns:
%    :
%
% Note:
%
% Example:
%
%    See also:

out=struct();
nvar=this.NumberOfVariables;
if nvar>1
    vnames=this.varnames;
    if ~isempty(vnames{1})
        for ivar=1:nvar
            S=struct('type','()','subs',{vnames(ivar)});
            obj=subsref(this,S);%<--obj=this(vnames{ivar});
            out.(vnames{ivar})=fanchart(obj,ci);
        end
        return
    end
end
if any(ci<0)
    error('confidence bands cannot be negative')
end
ci=sort(ci(:));
large=ci>=1;
ci(large)=ci(large)/100;
datax=squeeze(double(this));
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
out.date_numbers=this.date_numbers;