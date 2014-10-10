function [init,retcode]=kalman_initialization(T,R,steadystate,risk,transition_function,options)
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


% diffuse initialization for all elements in the state vector including
% stationary ones. This is what Waggoner and Zha do, but then they take
% a presample. The intuition, I guess, is that the filter eventually
% updates everything to the correct values. In some other cases, one
% may want set the presample to the number of unit roots as I have seen
% some place before... the drawback is that if the model has lots of
% unit roots and the sample is short...

if nargin==0
    if nargout>1
        error([mfilename,':: with no input argument, the number of output arguments cannot exceed 1'])
    end
    defaults=struct('kf_ergodic',true,... % formerly ergodic
        'kf_init_variance',[],... % formerly harvey_scale_factor
        'kf_presample',0,'kf_user_init',[]);
% % % % %         'kf_diffuse_min_presample',4,...
    % % %     'kf_diffuse_all',false,... % formerly diffuse_all
    lyap_options=lyapunov_equation();
    defaults=utils.miscellaneous.mergestructures(defaults,lyap_options);
    init=defaults;
    return
end
if nargin~=6
    error([mfilename,':: number of arguments must be 6'])
end

kf_ergodic=options.kf_ergodic;
kf_init_variance=options.kf_init_variance;
kf_presample=options.kf_presample;
kf_user_init=options.kf_user_init;

kf_diffuse_all=~isempty(kf_init_variance);

[init,retcode]=initialize_filter();

if ~retcode
    kf_presample=max(kf_presample,0);
    if kf_diffuse_all && kf_presample==0
        warning('Diffuse conditions detected with zero presample')
        warning('It is a good practive to have a presample period to initialize the filter under diffuse conditions')
    end
    init.start=kf_presample+1;
end

%--------------------------------------------------------------------------
    function [init,retcode]=initialize_filter()
        h=numel(T);
        endo_nbr=size(T{1},1);
        a0_given=numel(kf_user_init)>0;
        P0_given=numel(kf_user_init)>1;
        PAI00_given=numel(kf_user_init)>2;
        retcode=0;
        if PAI00_given
            PAI00=kf_user_init{3};
        else
            Q=transition_function(steadystate(:,1));
            [PAI00,retcode]=initial_markov_distribution(Q,kf_ergodic);
        end
        P0=[];
        a0=[];
        RR=cell(h,1);
        if ~retcode
            Tstar=0;
            Rstar=0;
            risk_star=zeros(endo_nbr,1);
            ss_star=0;
            for ireg=1:h
                Tstar=Tstar+PAI00(ireg)*T{ireg};
                Rstar=Rstar+PAI00(ireg)*R{ireg}(:,:,1);
                RR{ireg}=R{ireg}(:,:,1)*R{ireg}(:,:,1)';
                if ~isempty(risk)
                    risk_star=risk_star+PAI00(ireg)*risk{ireg};
                end
                ss_star=ss_star+PAI00(ireg)*steadystate{ireg};
            end
            if a0_given
                a0=kf_user_init{1};
            else
                a0=ss_star;
                if any(risk_star)
                    a0=a0+(eye(endo_nbr)-Tstar)\risk_star;
                end
            end
            if P0_given
                P0=kf_user_init{2};
            else
                if kf_diffuse_all
                    P0=kf_init_variance*eye(endo_nbr);
                else
                    [tmp,retcode]=lyapunov_equation(Tstar,Rstar*Rstar',options);
                    if ~retcode
                        P0=tmp;
                    end
                end
            end
            if ~retcode
                P0=repmat({P0},1,h);
                a0=repmat({a0},1,h);
            end
        end
        init=struct('a',{a0},'P',{P0},'PAI00',PAI00,'RR',{RR});
    end
end
