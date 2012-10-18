function plot_correlation_scatter(xparam,behave,parameter_names,pval_cutoff,r0,c0)

if ischar(parameter_names),parameter_names=cellstr(parameter_names);end

 parameter_names= parameter_names(:);

titel='Correlations of parameters in the Behavior sample';
scatter_correlations(xparam(:,behave),[],parameter_names,titel,pval_cutoff,r0,c0);

titel='Correlations of parameters in the Non-Behavior sample';
scatter_correlations(xparam(:,~behave),[],parameter_names,titel,pval_cutoff,r0,c0);

titel='Correlations of parameters in the Behavior(benchmark) and Non-Behavior samples';
scatter_correlations(xparam,behave,parameter_names,titel,pval_cutoff,r0,c0);

% npar=numel(parameter_names);
% 
% [Correlations,Pval]=corr(xparam(:,behave)');
% if npar<=15
%     disp('Correlations in the behavioral sample')
%     disp([cell(1,1),parameter_names;[parameter_names',num2cell(tril(Correlations))]])
% end
% 
% tmp=tril(Pval,-1);
% % detect significant correlations
% [rr,cc]=find(tmp<=pval_cutoff & tmp>0);
% NumberOfPlots=numel(rr);
% batch=[rr,cc];
% 
% disp(['Number of correlations seemingly significant= ',int2str(NumberOfPlots)])
% 
% titel='Correlations of parameters';
% nstar=r0*c0;
% nfig=ceil(NumberOfPlots/nstar);
% % SaveUnderName='ParameterCorrelations';
% for fig=1:nfig
%     if nfig>1
%         titelfig=[titel,' ',int2str(fig)];
%         %         SaveUnderName0=[SaveUnderName,'_',int2str(fig)];
%     else
%         %         SaveUnderName0=SaveUnderName;
%         titelfig=titel;
%     end
%     hfig=figure('name',titelfig);
%     [Remains,r,c]=number_of_rows_and_columns_in_figure(fig,NumberOfPlots,r0,c0);
%     for ii=1:min(nstar,Remains)
%         subplot(r,c,ii)
%         var_id=(fig-1)*nstar+ii;
%         p1=batch(var_id,1);
%         p2=batch(var_id,2);
%         d1=xparam(p1,~behave);
%         d2=xparam(p2,~behave);
%         d10=xparam(p1,behave);
%         d20=xparam(p2,behave);
%         plot(d1,d2,'.r','linewidth',2)
%         hold on
%         plot(d10,d20,'.b','linewidth',2)
%         axis tight
%         xlabel(parameter_names{p1},'interpreter','none')
%         ylabel(parameter_names{p2},'interpreter','none')
%         title([num2str(Correlations(p1,p2)),'(',num2str(Pval(p1,p2)),')'])
%         if ii==1
%             legend('Non-Behavior','behavior')
%         end
%     end
% end

function scatter_correlations(xparam,select,parameter_names,titel,pval_cutoff,r0,c0)
npar=numel(parameter_names);
if isempty(select)
    select=true(1,size(xparam,2));
end
[Correlations,Pval]=corr(xparam(:,select)');
if npar<=15
    disp('Correlations in the behavioral sample')
    disp([cell(1,1),parameter_names'
        [parameter_names,num2cell(tril(Correlations))]])
end
alternative=~all(select==true);

tmp=tril(Pval,-1);
% detect significant correlations
[rr,cc]=find(tmp<=pval_cutoff & tmp>0);
NumberOfPlots=numel(rr);
batch=[rr,cc];

disp(['Number of correlations seemingly significant= ',int2str(NumberOfPlots)])

nstar=r0*c0;
nfig=ceil(NumberOfPlots/nstar);
% SaveUnderName='ParameterCorrelations';
for fig=1:nfig
    if nfig>1
        titelfig=[titel,' ',int2str(fig)];
        %         SaveUnderName0=[SaveUnderName,'_',int2str(fig)];
    else
        %         SaveUnderName0=SaveUnderName;
        titelfig=titel;
    end
    hfig=figure('name',titelfig);
    [Remains,r,c]=number_of_rows_and_columns_in_figure(fig,NumberOfPlots,r0,c0);
    for ii=1:min(nstar,Remains)
        subplot(r,c,ii)
        var_id=(fig-1)*nstar+ii;
        p1=batch(var_id,1);
        p2=batch(var_id,2);
        d1=xparam(p1,select);
        d2=xparam(p2,select);
%         plot(d1,d2,'.r','linewidth',2)
        scatter(d1,d2,10,'b','filled','d')
        if alternative
            hold on
            d1=xparam(p1,~select);
            d2=xparam(p2,~select);
            scatter(d1,d2,10,'r','filled','d')
            if ii==1
                legend('Behavior','Non-Behavior')
            end
        end
        axis tight
        xlabel(parameter_names{p1},'interpreter','none')
        ylabel(parameter_names{p2},'interpreter','none')
        title([num2str(Correlations(p1,p2)),'(',num2str(Pval(p1,p2)),')'])
    end
end
