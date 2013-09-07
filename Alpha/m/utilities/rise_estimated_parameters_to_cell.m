function [Top,Middle,Bottom]=rise_estimated_parameters_to_cell(this)

nmod=numel(this);
estimated_parameters=[this.estimated_parameters];
npar=size(estimated_parameters,1);
no_item=' ';
Top=cell(6,1);
tabular=1;
Top{tabular}='\begin{tabular}{llrrr';
Top{2}='\hline\hline';
Top{3}=[no_item,'& \multicolumn{4}{c}{Prior}'];
Top{4}='\cline{2-5}';
Top{5}='Param & Distrib & Prob & low & high';
Top{6}='\hline';
if nmod>1
    for imod=1:nmod
        Top{tabular}=[Top{tabular},'rr']; % one for the skip and one for the results
        Top{3}=[Top{3},'&',no_item,'& Posterior'];
        index=int2str(5+2*imod);
        Top{4}=[Top{4},'\cline{',index,'-',index,'}'];
        Top{5}=[Top{5},'&',no_item,'& Mode(',int2str(imod),')'];
    end
else
    Top{tabular}=[Top{tabular},'rrrr']; % one for the skip, three for the results
    Top{3}=[Top{3},'&',no_item,'& \multicolumn{3}{c}{Posterior}'];
    Top{4}=[Top{4},'\cline{7-9}'];
    Top{5}=[Top{5},'&',no_item,'& Mode & 5\% & 95\%'];
end
Top{tabular}=[Top{tabular},'}'];
Top{5}=[Top{5},' \\'];
% now the parameters
Middle=cell(0,1);
for ipar=1:npar
    pname=estimated_parameters(ipar,1).name;
    pdistr=estimated_parameters(ipar,1).distribution;
    pprob=num2str(estimated_parameters(ipar,1).interval_probability);
    plow=num2str(estimated_parameters(ipar,1).plb);
    phigh=num2str(estimated_parameters(ipar,1).pub);
    newline=[pname,' & ',pdistr,' & ',pprob,' & ',plow,' & ',phigh];
    for imod=1:nmod
        pval=num2str(estimated_parameters(ipar,imod).mode);
        newline=[newline,' & ',no_item,' & ',pval]; %#ok<*AGROW>
        if nmod==1
            postlow=pval; % this is temporary
            posthigh=pval;
            % add the 5% and the 95% if available...
            newline=[newline,' & ',postlow,' & ',posthigh];
        end
    end
    newline=[newline,' \\ ']; % close the line
    Middle=[Middle;newline];
end
% closing
Bottom={'\hline';'\end{tabular}'};
