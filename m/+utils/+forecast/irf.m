function [irfs,retcode]=irf(y0,T,ss,states,which_shocks,Q,PAI,options)
% options=struct('simul_sig','simul_order','burn','nsimul','impulse','random','girf');

new_impulse=0;

[endo_nbr,nlags]=size(y0.y);

exo_nbr=numel(which_shocks);

nshocks=sum(which_shocks);

irfs=zeros(endo_nbr,nlags+options.nsteps,nshocks,options.nsimul);

iter=0;

det_vars=~which_shocks;

retcode=0;

for ishock=1:exo_nbr
    if det_vars(ishock)
        continue
    end
    iter=iter+1;
    shock_id=ishock;
    for isimul=1:options.nsimul
        if ~retcode
            [shocks]=utils.forecast.create_shocks(exo_nbr,shock_id,det_vars,options);
            if ~retcode
                [path1,~,retcode]=utils.forecast.multi_step(y0,ss,T,shocks,states,Q,PAI,options);
                if ~retcode
                    if options.girf
                        shocks=utils.forecast.replace_impulse(shocks,shock_id,options.k_future+1,new_impulse);
                        
                        [path2,~,retcode]=utils.forecast.multi_step(y0,ss,T,shocks,states,Q,PAI,options);
                        
                        path1=path1-path2;
                    end
                    
                    irfs(:,:,iter,isimul)=[y0.y,path1];
                end
            end
        end
    end
end