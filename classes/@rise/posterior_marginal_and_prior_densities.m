function posterior_marginal_and_prior_densities(obj)
correct_for_range=true;
simulation_folder=[obj.options.results_folder,filesep,'simulations'];
W = what(simulation_folder);
W=W.mat;
locs=find(strncmp('chain_',W,6));
if isempty(locs)
    error([mfilename,':: no simulations found'])
end
W=W(locs);
number_of_matrices=numel(W);
N=obj.options.discretize^2;
%=================
a=vertcat(obj.estimated_parameters.a);
b=vertcat(obj.estimated_parameters.b);
distr=vertcat(obj.estimated_parameters.distribution_code);
lb=vertcat(obj.estimated_parameters.lb);
ub=vertcat(obj.estimated_parameters.ub);
%=================
x0=vertcat(obj.estimated_parameters.mode);
f0=obj.log_post;
post_mode=x0;
f_post_mode=f0;
npar=size(x0,1);
r0=obj.options.graphics(1);
c0=obj.options.graphics(2);
nstar=r0*c0;
nfig=ceil(npar/nstar);
SaveUnderName0=[obj.options.results_folder,filesep,'graphs',filesep,'PriorsAndPosteriors'];
titel='priors and posterior marginal densities';
for fig=1:nfig
    [Remains,r,c]=number_of_rows_and_columns_in_figure(fig,npar,r0,c0);
    if nfig>1
        titelfig=[titel,' ',int2str(fig)];
		SaveUnderName=[SaveUnderName0,int2str(fig)];
    else
        titelfig=titel;
		SaveUnderName=SaveUnderName0;
    end
    hfig=figure('name',titelfig);
    for plt=1:min(nstar,Remains)
        par_id=(fig-1)*nstar+plt;
        all_vals=[];
        for m=1:number_of_matrices
            tmp=load([simulation_folder,filesep,W{m}]);
            Params=tmp.Params(par_id,:);
            if fig==1 && plt==1
                % try and locate the sampling posteriror mode
                fm=-tmp.minus_logpost_params;
                best=find(fm==max(fm),1,'first');
                if fm(best)>f_post_mode
                    post_mode=Params(:,best);
                    f_post_mode=fm(best);
                end
            end
            all_vals=[all_vals;Params(:)]; %#ok<AGROW>
        end
        
        figure(hfig)
		mm=mean(all_vals);
        xmin = min(all_vals);
        xmax = max(all_vals);
%         grid = transpose(linspace(xmin,xmax,N));
% 		[F,XI]=ksdensity(all_vals,grid);
		[F,XI]=distributions.kernel_density(all_vals,[],[],'normal',N);
        [x_mode,x_mode_id]=find_nearest(XI,x0(par_id));
        [x_post_mode,x_post_mode_id]=find_nearest(XI,post_mode(par_id));
        [x_mm,x_mm_id]=find_nearest(XI,mm);
        x_prior=linspace(lb(par_id),ub(par_id),N);
        x_prior=x_prior(:);
        f_prior=prior_density(distr(par_id),x_prior,a(par_id),b(par_id));
        if correct_for_range
            % give it the same range as F
            if max(f_prior)==min(f_prior)
                ratio=.5;
            else
                ratio=(f_prior-min(f_prior))/(max(f_prior)-min(f_prior));
            end
            f_prior=min(F)+ratio*(max(F)-min(F));
        end
        subplot(r,c,plt)
        plot(XI,F,'LineStyle','-','Color','b',...
            'LineWidth',2.5), hold on
        plot(x_prior,f_prior,'LineStyle','-','Color','green',...
            'LineWidth',2.5), hold on
        plot([x_mm x_mm], [0 F(x_mm_id)],'LineStyle',':',...
            'Color','black','LineWidth',2.5 ),
        plot([x_mode x_mode], [0 F(x_mode_id)],'LineStyle',':',...
            'Color','green','LineWidth',2.5 ),
        plot([x_post_mode x_post_mode], [0 F(x_post_mode_id)],'LineStyle',':',...
            'Color','red','LineWidth',2.5 ),
        hold off
        xlow=min(xmin,min(x_prior));
        xhigh=max(xmax,max(x_prior));
        axis([xlow xhigh 0 1.4*max(F)]);
        title(obj.estimated_parameters(par_id).tex_name,...
            'FontSize',12,'FontWeight','bold')%'interpreter','none',
        if plt==1
            legend('post density','prior density','mean','mode','post_mode')
        end
    end
	saveas(hfig,[SaveUnderName,'.pdf'])
    saveas(hfig,[SaveUnderName,'.fig'])
    saveas(hfig,[SaveUnderName,'.eps'])
end
