function [irfs,retcode]=irf(y0,T,ss,state_vars_location,which_shocks,det_vars,options)
% H1 line
%
% ::
%
%
% Args:
%
% Returns:
%    :
%
% Note:
%
% Example:
%
%    See also:

new_impulse=0;

[endo_nbr,nlags]=size(y0.y);

exo_nbr=numel(which_shocks);

nshocks=sum(which_shocks);

if any(which_shocks & det_vars)
    
    error('trying to compute irfs for deterministic shocks')
    
end

irfs=zeros(endo_nbr,nlags+options.nsteps,nshocks,options.nsimul);

iter=0;

retcode=0;
% orig_shocks=y0.econd.data(:,:,1);
for ishock=1:exo_nbr
    
    if ~which_shocks(ishock)
        
        continue
        
    end
    
    iter=iter+1;
    
    shock_id=ishock;
    
    for isimul=1:options.nsimul
        
        if ~retcode
            
            shocks=utils.forecast.create_shocks(exo_nbr,shock_id,det_vars,options);
            % make sure that impulses that are inherited stay on in the
            % reference simulation. This may imply over-riding the impulse
            % itself if it happens to be on the path of an inherited shock
%             %--------------------------------------------------------------
%             inherited_shocks=orig_shocks~=0;
%             shocks(inherited_shocks)=orig_shocks(inherited_shocks);
            y0.econd.data=shocks(:,:,ones(3,1));
            
            if ~retcode
                
                [sim1,states1,retcode]=utils.forecast.multi_step(y0,ss,T,state_vars_location,options);
                
                if ~retcode
                    
                    path1=[y0.y,sim1];
                    
                    if options.girf
                        
                        if isimul==1
                            
                            path2=path1;
                            
                        end
                        
                        y02=y0;
                        % set to 0(new_impulse) the location
                        % corresponding to the impulse
                        %------------------------------------------
                        shocks=utils.forecast.replace_impulse(shocks,shock_id,...
                            options.k_future+1,new_impulse);
                        % set to 0 the location corresponding to the
                        % inherited shocks
                        %------------------------------------------
%                         shocks(inherited_shocks)=0;
                        y02.econd.data=shocks(:,:,ones(3,1));
                        
                        states2=states1;
                        
                        if options.girf_regime_uncertainty
                            
                            states2=nan(size(states1));
                            
                        end
                        
                        y02.rcond.data=states2;
                        % ensure that the shocks are not updated in the
                        % alternative scenario
                        %-----------------------------------------------
                        [sim2,~,retcode]=utils.forecast.multi_step(y02,ss,T,state_vars_location,options);
                        
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