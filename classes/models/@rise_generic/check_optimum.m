function [db,fighandles]=check_optimum(obj,plotit,varlist)%,r0,c0,N
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
% More About
% ------------
%
% Examples
% ---------
%
% See also: 

if isempty(obj)
    if nargout>1
        error([mfilename,':: when the object is emtpy, nargout must be at most 1'])
    end
    db=struct();
    return
end
if nargin<3
    varlist=[];
    if nargin<2
        plotit=[];
    end
end
AllNames={obj.estimation.priors.name};
if isempty(varlist)
    varlist=AllNames;
end
if isempty(plotit)
    plotit=true;
end

if ischar(varlist)
    varlist=cellstr(varlist);
end
list_locs=locate_variables(varlist,AllNames);

r0=obj.options.graphics(1);
c0=obj.options.graphics(2);
N=obj.options.prior_discretize;

LB=vertcat(obj.estimation.priors.lower_bound);
UB=vertcat(obj.estimation.priors.upper_bound);
xmode=obj.estimation.posterior_maximization.mode;
varnames={obj.estimation.priors.tex_name};
interpreter='none';
if strcmp(varnames{1}(1),'\')
    interpreter='latex';
end

nvar=numel(varlist);
nstar=r0*c0;
nfig=ceil(nvar/nstar);
fighandles=nan(nfig,1);
obj.options.kf_filtering_level=0; % do not filter

disp([mfilename,':: using serial code'])
npar=numel(list_locs);
db=cell(2,npar);
mainfunc=@do_one_parameter;
varnames={obj.estimation.priors(list_locs).name};
parfor ipar=1:npar
    var_id=list_locs(ipar);
    vname=varnames{var_id}; %#ok<PFBNS>
    db(ipar,:)={vname,mainfunc(var_id)};
end

fields=db(:,1);
db=cell2struct(db(:,2),fields,1);

if plotit
    titel='Check plots';
    for fig=1:nfig
        if nfig>1
            titelfig=[titel,' ',int2str(fig)];
        else
            titelfig=titel;
        end
        fighandles(fig)=figure('name',titelfig);
        [Remains,r,c]=number_of_rows_and_columns_in_figure(fig,nvar,r0,c0);
        for ii=1:min(nstar,Remains)
            subplot(r,c,ii)
            var_id=list_locs((fig-1)*nstar+ii);
            vname=obj.estimation.priors(var_id).name;
            do_one_plot(db.(vname))
        end
    end
end

    function pp=do_one_parameter(var_id)
        pp=struct();
        vtexname=obj.estimation.priors(var_id).tex_name;
        pp.tex_name=vtexname;
        pp.mode=xmode(var_id);
        pp.log_post_mode=obj.estimation.posterior_maximization.log_post;
        pp.log_lik_mode=obj.estimation.posterior_maximization.log_lik;
        low = max(LB(var_id),0.8*xmode(var_id));
        high = min(UB(var_id),1.2*xmode(var_id));
        pp.x=sort([linspace(low,high,N),xmode(var_id)]);
        posj=find(abs(pp.x-xmode(var_id))==min(abs(pp.x-xmode(var_id))),1,'first');
        pp.log_post=zeros(1,N+1);
        pp.log_lik=zeros(1,N+1);
        for jj=1:N+1
            if jj~=posj
                pj=xmode;
                pj(var_id)=pp.x(jj);
                [pp.log_post(jj),pp.log_lik(jj)]=log_posterior_kernel(obj,pj);
            else
                pp.log_post(jj)=pp.log_post_mode;
                pp.log_lik(jj)=pp.log_lik_mode;
            end
        end
        pp.log_post(pp.log_lik<=-obj.options.estim_penalty)=nan;
        pp.log_lik(pp.log_lik<=-obj.options.estim_penalty)=nan;
    end

    function do_one_plot(pp)
        low_f=min(min([pp.log_post,pp.log_lik]));
        high_f=max(max([pp.log_post,pp.log_lik]));
        posj=find(abs(pp.x-pp.mode)==min(abs(pp.x-pp.mode)),1,'first');
        plot(pp.x,pp.log_post,...
            pp.x,pp.log_lik,...
            [pp.x(posj),pp.x(posj)],[low_f,high_f],...
            'linewidth',1.5)
        title(pp.tex_name,'interpreter',interpreter) %
        if ii==1
            legend({'log post','log lik','mode'},'location','SW','orientation','horizontal')
        end
        if any(isnan(pp.log_post))
            hold on
            locs=find(isnan(pp.log_post));
            plot(pp.x(locs)',min(pp.log_post)*ones(numel(locs),1),'.r','markersize',10)
            hold off
        end
        axis tight % xlim([min(pp.x),max(pp.x)])
    end
end

% %     SaveUnderName0=[obj.folders_paths.graphs,filesep,'CheckPlots'];
% %             SaveUnderName=[SaveUnderName0,int2str(fig)];
% % %             SaveUnderName=SaveUnderName0;
        %     saveas(hfig,[SaveUnderName,'.pdf'])
        %     saveas(hfig,[SaveUnderName,'.fig'])
        %     saveas(hfig,[SaveUnderName,'.eps'])
