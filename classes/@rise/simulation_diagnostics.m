function simulation_diagnostics(obj)
simulation_folder=[obj.options.results_folder,filesep,'simulations'];
W = what(simulation_folder);
W=W.mat;
locs=find(strncmp('chain_',W,6));
if isempty(locs)
    error([mfilename,':: no simulations found'])
end
W=W(locs);
% first determine the number of chains and the number  of matrices in each
% chain.
number_of_parallel_chains=0;
number_of_matrices=0;
for ii=1:numel(W)
    chain_name=W{ii};
    underscores=strfind(chain_name,'_');
    dot_id=strfind(chain_name,'.');
    chain_nbr=str2double(chain_name(underscores(1)+1:underscores(2)-1));
    number_of_parallel_chains=max(number_of_parallel_chains,chain_nbr);
    matrix_nbr=str2double(chain_name(underscores(2)+1:dot_id-1));
    number_of_matrices=max(number_of_matrices,matrix_nbr);
end

npar=size(obj.estimated_parameters,1);
r0=obj.options.graphics(1);
c0=obj.options.graphics(2);
nstar=r0*c0;
nfig=ceil(npar/nstar);
test_names={'RecursiveMeans','RecursiveVariances','PSRF'};
titel={'recursive means','recursive variances','potential scale reduction'};
tmp=[obj.options.results_folder,filesep,'graphs',filesep];
SaveUnderName0=cell(1,3);
titelfig=cell(1,3);
for ii=1:3
	SaveUnderName0{ii}=[tmp,test_names{ii}];
end
SaveUnderName=cell(1,3);
hfig=cell(1,3);
for fig=1:nfig
    [Remains,r,c]=number_of_rows_and_columns_in_figure(fig,npar,r0,c0);
	for ii=1:3
	    hfig{ii}=figure('name',[titel{ii},int2str(fig)]);
	    if nfig>1
	        titelfig{ii}=[titel{ii},' ',int2str(fig)];
			SaveUnderName{ii}=[SaveUnderName0{ii},int2str(fig)];
	    else
	        titelfig{ii}=titel{ii};
			SaveUnderName{ii}=SaveUnderName0{ii};
	    end
	end
    for plt=1:min(nstar,Remains)
        par_id=(fig-1)*nstar+plt;
        par_mean=cell(number_of_parallel_chains,1);
        par_var=cell(number_of_parallel_chains,1);
        all_vals=cell(number_of_parallel_chains,1);
        for pc=1:number_of_parallel_chains
            iter=0;
            m0=0;
            V0=0;
            for m=1:number_of_matrices
                tmp=load([simulation_folder,filesep,'chain_',int2str(pc),'_',int2str(m)]);
                Params=tmp.Params(par_id,:);
                nvals=numel(Params);
                all_vals{pc}=[all_vals{pc},Params];
                for ii=1:nvals
                    iter=iter+1;
                    [m0,V0]=recursive_moments(m0,V0,Params(ii),iter);
                    par_mean{pc}=[par_mean{pc},m0];
                    par_var{pc}=[par_var{pc},V0];
                end
            end
        end
        par_mean=cell2mat(par_mean);
        par_var=cell2mat(par_var);
        all_vals=cell2mat(all_vals);
        
        figure(hfig{1})
        subplot(r,c,plt)
        plot(par_mean','linewidth',1.5)
        title(obj.estimated_parameters(par_id).name,'interpreter','none')
        
        figure(hfig{2})
        subplot(r,c,plt)
        plot(par_var','linewidth',1.5)
        title(obj.estimated_parameters(par_id).name,'interpreter','none')
        
        % potential scale reduction
        R=potential_scale_reduction(all_vals);
        figure(hfig{3})
        subplot(r,c,plt)
        plot(R,'linewidth',1.5)
        title(obj.estimated_parameters(par_id).name,'interpreter','none')
    end
    figure(hfig{1}),sup_label('recursive means','y');
    figure(hfig{2}),sup_label('recursive variances','y');
    figure(hfig{3}),sup_label('potential scale reduction','y');
	for ii=1:3
		saveas(hfig{ii},[SaveUnderName{ii},'.pdf'])
	    saveas(hfig{ii},[SaveUnderName{ii},'.fig'])
	    saveas(hfig{ii},[SaveUnderName{ii},'.eps'])
	end
end
