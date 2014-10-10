function retcode=prior_plots(obj)
% H1 line
%
% Syntax
% -------
% ::
%
% Inputs
% -------
%
% Outputs
% --------
%
% Description
% ------------
%
% Examples
% ---------
%
% See also: 


if isempty(obj)
	% For the computation of check plots, priors and posteriors
    retcode=struct('prior_discretize',20);
    return
end
if nargout
    retcode=0;
end
default_cutoff=1e-4;
nobj=numel(obj);

r0=obj(1).options.graphics(1);
c0=obj(1).options.graphics(2);
N=40*obj(1).options.prior_discretize;
if isempty(obj(1).estimation)
    error([mfilename,':: No estimated parameters information provided'])
end

SaveUnderName0=cell(1,nobj);%
npriors=numel(obj(1).estimation.priors);
a=nan(npriors,nobj);
b=nan(npriors,nobj);
% distr=cell(npriors,nobj);
lb=nan(npriors,nobj);
ub=nan(npriors,nobj);
for ii=1:nobj
    a(:,ii)=[obj(ii).estimation.priors.a]';
    b(:,ii)=[obj(ii).estimation.priors.b]';
%     distr(:,ii)=transpose({obj(ii).estimated_parameters.distribution});
    lb(:,ii)=[obj(ii).estimation.priors.lower_bound]';
    ub(:,ii)=[obj(ii).estimation.priors.upper_bound]';
    SaveUnderName0{ii}=[obj(ii).options.results_folder,filesep,'graphs',filesep,'Priors'];
end

nvar=numel(obj(1).estimation.priors);
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
            ldens=distributions.(obj(jj).estimation.priors(var_id).prior_distrib)();
            aa=obj(jj).estimation.priors(var_id).a;
            bb=obj(jj).estimation.priors(var_id).b;
            cc=[];%obj(jj).estimation.priors(var_id).c;
            dd=[];%obj(jj).estimation.priors(var_id).d;
            f(:,jj)=exp(...
                ldens(x(:,jj),aa,bb,cc,dd)...
                );
%             %=========================
%             significant=f(:,jj)>default_cutoff;
%             while sum(significant)<N/3
%                 keyboard
%                 x(:,jj)=0.1*x(:,jj);
%                 f(:,jj)=exp(ldens(x(:,jj),aa,bb,cc,dd));
%                 significant=f(:,jj)>default_cutoff;
%             end
%             %=========================
        end
        subplot(r,c,ii)
        plot(x,f)
        title(obj(1).estimation.priors(var_id).tex_name,'Interp','tex')% ,'interpreter','none'
        low_f=utils.stat.nanmin(f(:));
        high_f=utils.stat.nanmax(f(:));
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