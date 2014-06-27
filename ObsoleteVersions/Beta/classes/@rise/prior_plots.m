function prior_plots(obj)

default_cutoff=1e-4;
nobj=numel(obj);

r0=obj(1).options.graphics(1);
c0=obj(1).options.graphics(2);
N=40*obj(1).options.discretize;
if isempty(obj(1).estimated_parameters)
    error([mfilename,':: No estimated parameters information provided'])
end

SaveUnderName0=cell(1,nobj);%
a=nan(numel(obj(1).estimated_parameters),nobj);
b=nan(numel(obj(1).estimated_parameters),nobj);
% distr=cell(numel(obj(1).estimated_parameters),nobj);
lb=nan(numel(obj(1).estimated_parameters),nobj);
ub=nan(numel(obj(1).estimated_parameters),nobj);
for ii=1:nobj
    a(:,ii)=vertcat(obj(ii).estimated_parameters.a);
    b(:,ii)=vertcat(obj(ii).estimated_parameters.b);
%     distr(:,ii)=transpose({obj(ii).estimated_parameters.distribution});
    lb(:,ii)=vertcat(obj(ii).estimated_parameters.lb);
    ub(:,ii)=vertcat(obj(ii).estimated_parameters.ub);
    SaveUnderName0{ii}=[obj(ii).options.results_folder,filesep,'graphs',filesep,'Priors'];
end

nvar=size(obj(1).estimated_parameters,1);
nstar=r0*c0;
nfig=ceil(nvar/nstar);

titel='Prior plots';
for fig=1:nfig
    if nfig>1
        titelfig=[titel,' ',int2str(fig)];
    else
        titelfig=titel;
    end
    hfig=figure('name',titelfig); % I will need a handle if I ever want to save the figure.
    [Remains,r,c]=number_of_rows_and_columns_in_figure(fig,nvar,r0,c0);
    for ii=1:min(nstar,Remains)
        var_id=(fig-1)*nstar+ii;
        x=nan(N,nobj);
        f=nan(N,nobj);
        for jj=1:nobj
            x(:,jj)=linspace(lb(var_id,jj),ub(var_id,jj),N);
            ldens=distributions.(obj(jj).estimated_parameters(var_id).distribution)();
            aa=obj(jj).estimated_parameters(var_id).a;
            bb=obj(jj).estimated_parameters(var_id).b;
            cc=obj(jj).estimated_parameters(var_id).c;
            dd=obj(jj).estimated_parameters(var_id).d;
            f(:,jj)=exp(...
                ldens(x(:,jj),aa,bb,cc,dd)...
                );
            %=========================
            significant=f(:,jj)>default_cutoff;
            while sum(significant)<N/3
                x(:,jj)=0.1*x(:,jj);
                f(:,jj)=exp(ldens(x(:,jj),aa,bb,cc,dd));
                significant=f(:,jj)>default_cutoff;
            end
            %=========================
        end
        subplot(r,c,ii)
        plot(x,f)
        title(obj(1).estimated_parameters(var_id).tex_name,'Interp','tex')% ,'interpreter','none'
        low_f=nanmin(nanmin(f));
        high_f=nanmax(nanmax(f));
        axis([min(min(x)),max(max(x)),low_f-abs(low_f)/1000,high_f+abs(high_f)/1000])% tight
        if nobj>1 && ii==1
            legend({obj.filename})
        end
    end
    for jj=1:nobj
        SaveUnderName=SaveUnderName0{jj};
        if nfig>1
            SaveUnderName=[SaveUnderName,int2str(fig)]; %#ok<AGROW>
        end
        saveas(hfig,[SaveUnderName,'.pdf'])
        saveas(hfig,[SaveUnderName,'.fig'])
        saveas(hfig,[SaveUnderName,'.eps'])
    end
end