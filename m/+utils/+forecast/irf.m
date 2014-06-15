function [irfs,retcode]=irf(y0,T,ss,which_shocks,options)

new_impulse=0;

[endo_nbr,nlags]=size(y0.y);

exo_nbr=numel(which_shocks);

nshocks=sum(which_shocks);

irfs=zeros(endo_nbr,nlags+options.nsteps,nshocks,options.nsimul);

iter=0;

det_vars=~which_shocks;

retcode=0;
path1=nan(endo_nbr,nlags+options.nsteps);
path1(:,1:nlags)=y0.y;
for ishock=1:exo_nbr
    if det_vars(ishock)
        continue
    end
    iter=iter+1;
    shock_id=ishock;
    for isimul=1:options.nsimul
        if ~retcode
            options.shocks=utils.forecast.create_shocks(exo_nbr,shock_id,det_vars,options);
            if ~retcode
                [sim1,~,retcode]=utils.forecast.multi_step(y0,ss,T,options);
                if ~retcode
                    path1(:,nlags+1:end)=sim1;
                    if options.girf
                        if isimul==1
                            path2=path1;
                        end
                        options.shocks=utils.forecast.replace_impulse(options.shocks,shock_id,options.k_future+1,new_impulse);
                        
                        [sim2,~,retcode]=utils.forecast.multi_step(y0,ss,T,options);
                        if ~retcode
                            path2(:,nlags+1:end)=sim2;
                            path1=path1-path2;
                        end
                    end
                    irfs(:,:,iter,isimul)=path1;
                end
            end
        end
    end
end