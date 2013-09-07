function check_optimum(obj,varlist)%,r0,c0,N
if nargin==1
    varlist=[];
end
AllNames={obj.estimated_parameters.name};
if isempty(varlist)
    varlist=AllNames;
end
if ischar(varlist)
    varlist=cellstr(varlist);
end
list_locs=locate_variables(varlist,AllNames);
estimated_parameters=obj.estimated_parameters;

r0=obj.options.graphics(1);
c0=obj.options.graphics(2);
N=obj.options.discretize;

LB=vertcat(estimated_parameters.lb);
UB=vertcat(estimated_parameters.ub);
xmode=vertcat(estimated_parameters.mode);
varnames={estimated_parameters.tex_name};
interpreter='none';
if strcmp(varnames{1}(1),'\')
    interpreter='latex';
end

nvar=numel(varlist);
nstar=r0*c0;
nfig=ceil(nvar/nstar);

titel='Check plots';
SaveUnderName0=[obj.options.results_folder,filesep,'graphs',filesep,'CheckPlots'];

% if exist('matlabpool.m','file') && matlabpool('size')>0
%     disp([mfilename,':: using parallel code: the graphs will not be shown'])
%     parfor fig=1:nfig
%         if nfig>1
%             titelfig=[titel,' ',int2str(fig)];
%             SaveUnderName=[SaveUnderName0,int2str(fig)];
%         else
%             titelfig=titel;
%             SaveUnderName=SaveUnderName0;
%         end
%         hfig=figure('name',titelfig); % I will need a handle if I ever want to save the figure.
%         [Remains,r,c]=number_of_rows_and_columns_in_figure(fig,nvar,r0,c0);
%         for ii=1:min(nstar,Remains)
%             var_id=list_locs((fig-1)*nstar+ii);
%             low = max(LB(var_id),0.8*xmode(var_id)); %#ok<*PFBNS>
%             high = min(UB(var_id),1.2*xmode(var_id));
%             xj=sort([linspace(low,high,N),xmode(var_id)]);
%             posj=find(abs(xj-xmode(var_id))==min(abs(xj-xmode(var_id))),1,'first');
%             log_post=zeros(N+1,1);
%             log_lik=zeros(N+1,1);
%             for jj=1:N+1
%                 if jj~=posj
%                     pj=xmode;
%                     pj(var_id)=xj(jj);
%                     [log_post(jj),log_lik(jj)]=log_posterior_kernel(obj,pj,false); % do not filter
%                 else
%                     log_post(jj)=obj.log_post;
%                     log_lik(jj)=obj.log_lik;
%                 end
%             end
%             log_post(log_lik<=-obj.options.Penalty)=nan;
%             log_lik(log_lik<=-obj.options.Penalty)=nan;
%             subplot(r,c,ii)
%             plottitle=deblank(varnames{var_id});
%             low_f=min(min([log_post,log_lik]));
%             high_f=max(max([log_post,log_lik]));
%             plot(xj',log_post,...
%                 xj',log_lik,...
%                 [xj(posj),xj(posj)],[low_f,high_f],...
%                 'linewidth',1.5)
%             title(plottitle,'interpreter',interpreter) %
%             if ii==1
%                 legend('log post','log lik','mode')
%             end
%             if any(isnan(log_post))
%                 hold on
%                 locs=find(isnan(log_post));
%                 plot(xj(locs)',min(log_post)*ones(numel(locs),1),'.r','markersize',10)
%                 hold off
%             end
%             %         axis([min(xj),max(xj),low_f-abs(low_f)/1000,high_f+abs(high_f)/1000])% tight
%             xlim([min(xj),max(xj)])
%         end
%         saveas(hfig,[SaveUnderName,'.pdf'])
%         saveas(hfig,[SaveUnderName,'.fig'])
%         saveas(hfig,[SaveUnderName,'.eps'])
%     end
% else
disp([mfilename,':: using serial code'])
for fig=1:nfig
    if nfig>1
        titelfig=[titel,' ',int2str(fig)];
        SaveUnderName=[SaveUnderName0,int2str(fig)];
    else
        titelfig=titel;
        SaveUnderName=SaveUnderName0;
    end
    hfig=figure('name',titelfig); % I will need a handle if I ever want to save the figure.
    [Remains,r,c]=number_of_rows_and_columns_in_figure(fig,nvar,r0,c0);
    for ii=1:min(nstar,Remains)
        var_id=list_locs((fig-1)*nstar+ii);
        low = max(LB(var_id),0.8*xmode(var_id));
        high = min(UB(var_id),1.2*xmode(var_id));
        xj=sort([linspace(low,high,N),xmode(var_id)]);
        posj=find(abs(xj-xmode(var_id))==min(abs(xj-xmode(var_id))),1,'first');
        log_post=zeros(N+1,1);
        log_lik=zeros(N+1,1);
        for jj=1:N+1
            if jj~=posj
                pj=xmode;
                pj(var_id)=xj(jj);
                [log_post(jj),log_lik(jj)]=log_posterior_kernel(obj,pj,false); % do not filter
            else
                log_post(jj)=obj.log_post;
                log_lik(jj)=obj.log_lik;
            end
        end
        log_post(log_lik<=-obj.options.Penalty)=nan;
        log_lik(log_lik<=-obj.options.Penalty)=nan;
        subplot(r,c,ii)
        plottitle=deblank(varnames{var_id});
        low_f=min(min([log_post,log_lik]));
        high_f=max(max([log_post,log_lik]));
        plot(xj',log_post,...
            xj',log_lik,...
            [xj(posj),xj(posj)],[low_f,high_f],...
            'linewidth',1.5)
        title(plottitle,'interpreter',interpreter) %
        if ii==1
            legend('log post','log lik','mode')
        end
        if any(isnan(log_post))
            hold on
            locs=find(isnan(log_post));
            plot(xj(locs)',min(log_post)*ones(numel(locs),1),'.r','markersize',10)
            hold off
        end
        %         axis([min(xj),max(xj),low_f-abs(low_f)/1000,high_f+abs(high_f)/1000])% tight
        xlim([min(xj),max(xj)])
    end
    saveas(hfig,[SaveUnderName,'.pdf'])
    saveas(hfig,[SaveUnderName,'.fig'])
    saveas(hfig,[SaveUnderName,'.eps'])
end
% end

