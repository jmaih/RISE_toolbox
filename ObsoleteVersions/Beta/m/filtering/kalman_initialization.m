function [a0,P0,PAI00,start,retcode]=kalman_initialization(T,R,transition_matrix,options)
% There is no documentation of this function yet.

% diffuse initialization for all elements in the state vector including
% stationary ones. This is what Waggoner and Zha do, but then they take
% a presample. The intuition, I guess, is that the filter eventually
% updates everything to the correct values. In some other cases, one
% may want set the presample to the number of unit roots as I have seen
% some place before... the drawback is that if the model has lots of
% unit roots and the sample is short...

defaults=struct('kf_ergodic',true,... % formerly ergodic
    'kf_init_variance',[],... % formerly harvey_scale_factor
    'kf_diffuse_min_presample',4,...
    'kf_presample',0);
% % %     'kf_diffuse_all',false,... % formerly diffuse_all
lyap_options=lyapunov_equation();
defaults=mergestructures(defaults,lyap_options);
if nargin==0
    if nargout>1
        error([mfilename,':: with no input argument, the number of output arguments cannot exceed 1'])
    end
    a0=defaults;
    return
end
if nargin<4
    options=[];
    if nargin<3
        error([mfilename,':: this function requires at least 3 input arguments'])
    end
elseif nargin>7
    error([mfilename,':: number of arguments cannot exceed 4'])
end

init_fields=fieldnames(defaults);
for ii=1:numel(init_fields)
    v=init_fields{ii};
    if isfield(options,v)
        defaults.(v)=options.(v);
    end
end
kf_ergodic=defaults.kf_ergodic;
kf_init_variance=defaults.kf_init_variance;
kf_diffuse_min_presample=defaults.kf_diffuse_min_presample;
kf_presample=defaults.kf_presample;

kf_diffuse_all=~isempty(kf_init_variance);

[a0,P0,PAI00,retcode]=initialize_filter();

start=[]; 
if ~retcode
    kf_presample=max(kf_presample,0);
    if kf_diffuse_all
        kf_presample=max(kf_presample,kf_diffuse_min_presample);
    end
    start=kf_presample+1;
end

%--------------------------------------------------------------------------
    function [a0,P0,PAI00,retcode]=initialize_filter()
        [endo_nbr,~,h]=size(T);
        [PAI00,retcode]=initial_markov_distribution(transition_matrix,kf_ergodic);
        P0=[];
        a0=[];
        if ~retcode
            if kf_diffuse_all
                P0=kf_init_variance*eye(endo_nbr);
            else
                Tstar=0;
                Rstar=0;
                for ireg=1:h
                    Tstar=Tstar+PAI00(ireg)*T(:,:,ireg);
                    Rstar=Rstar+PAI00(ireg)*R(:,:,ireg);
                end
                [tmp,retcode]=lyapunov_equation(Tstar,Rstar*Rstar',options);
                if ~retcode
                    P0=tmp;
                end
            end
            if ~retcode
                P0=P0(:,:,ones(1,h),ones(1,h));
                a0=zeros(endo_nbr,h,h);
            end
        end
    end
end
