function [irfs,retcode]=irf(y0,T,ss,state_vars_location,which_shocks,det_vars,options)
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
orig_options=options;
for ishock=1:exo_nbr
    if ~which_shocks(ishock)
        continue
    end
    iter=iter+1;
    shock_id=ishock;
    for isimul=1:options.nsimul
        if ~retcode
            options.shocks=utils.forecast.create_shocks(exo_nbr,shock_id,det_vars,orig_options);
            % make sure that impulses that are inherited stay on in the
            % reference simulation. This may imply over-riding the impulse
            % itself if it happens to be on the path of an inherited shock
            %--------------------------------------------------------------
            inherited_shocks=orig_options.shocks~=0;
            options.shocks(inherited_shocks)=orig_options.shocks(inherited_shocks);
            if ~retcode
                [sim1,states1,retcode]=utils.forecast.multi_step(y0,ss,T,state_vars_location,options);
                if ~retcode
                    path1=[y0.y,sim1];
                    if options.girf
                        if isimul==1
                            path2=path1;
                        end
                        options2=options;
                        % set to 0(new_impulse) the location
                        % corresponding to the impulse
                        %------------------------------------------
                        options2.shocks=utils.forecast.replace_impulse(options.shocks,shock_id,options.k_future+1,new_impulse);
                        % set to 0 the location corresponding to the
                        % inherited shocks
                        %------------------------------------------
                        options2.shocks(inherited_shocks)=0;
                        options2.states=states1;
                        % ensure that the shocks are not updated in the
                        % alternative scenario
                        %-----------------------------------------------
                        options2.simul_do_update_shocks=false;
                        [sim2,~,retcode]=utils.forecast.multi_step(y0,ss,T,state_vars_location,options2);
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