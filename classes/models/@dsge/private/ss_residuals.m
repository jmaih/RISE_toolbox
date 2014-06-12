function [r,Jac,retcode]=ss_residuals(ss_0,resid_func,func_jac,obj)

done=obj.is_optimal_policy_model || obj.is_imposed_steady_state;
retcode=0;
Jac=[];
r=0*ss_0(:);

if ~done
    % input list is always 'y'  'x'  'ss'  'param'  'sparam' 'def'  's0'  's1'
    x_ss=zeros(sum(obj.exogenous.number),1);
    def=obj.solution.definitions;
    pp=obj.parameter_values;
    number_of_regimes=obj.markov_chains.regimes_number;
    ss_0=reshape(ss_0,[],number_of_regimes);
    [TransMat,retcode]=compute_steady_state_transition_matrix(obj,ss_0(:,1),pp);
    if ~retcode
        % compute the unique steady state based on the ergodic distribution
        %------------------------------------------------------------------
        if obj.is_unique_steady_state
            [pp_i,def_i,retcode]=ergodic_parameters(TransMat.Qinit,def,pp);
        end
        
        % compute the residuals
        %----------------------
        r=ss_0;
        % check that the steady states computed above are actually the steady
        % states. If not, then compute the steady state(s)
        for ii=1:number_of_regimes
            % if the initial guess solves the steady state then proceed. Else
            % try and improve the initial guess through fsolve.
            if ~obj.is_unique_steady_state
                def_i=def{ii};
                pp_i=pp(:,ii);
            end
            y_=ss_0(:,ii);
            ss_=y_;
            % input list is always 'y'    'x'    'ss'    'param'    'sparam'    'def'    's0'    's1'
            r(:,ii)=utils.code.evaluate_functions(resid_func,y_,x_ss,ss_,pp_i,[],def_i,[],[]);
            if nargout>1
                Jac=func_jac(y_,x_ss,ss_,pp_i,[],def_i,[],[]);
            end
            if  obj.is_unique_steady_state && ...
                    ii==1 && number_of_regimes>1
                test=bsxfun(@minus,ss_0(:,1),ss_0(:,2:end));
                if all(abs(test(:))<1e-7)
                    r=r(:,ii*ones(1,number_of_regimes));
                    break
                end
            end
        end
        r=r(:);
    end
end
if retcode && obj(1).options.debug
    utils.error.decipher(retcode)
end
end