function [ld,user_info,retcode]=ctw_endogenous_priors(obj,user_info,targets,prior_sample)

if nargin<4
    prior_sample=[];
    if nargin<3
        targets=[];
    end
end

if isempty(user_info)
    initialize();
end

[P,retcode]=theoretical_autocovariances(obj,0);

if retcode
    ld=nan;
else
    ld=log_density();
end

    function ld=log_density() 
        Ftheta=diag(P(user_info.ids,user_info.ids));
        FF=user_info.Fhat-Ftheta;
        ld=.5*user_info.ny*log(user_info.T/(2*pi))...
            -.5*log(det(user_info.Shat))...
            -.5*user_info.T*(FF'/user_info.Shat)*FF;
    end

    function initialize()
        user_info=struct();
        user_info.ll=2;
        user_info.theta=1;
        if isempty(targets)
            targets=obj.observables.name;
        elseif ischar(targets)
            targets=cellstr(targets);
        end
        user_info.targets=targets;
        %-----------------------------------------------
        data=obj.data.y(:,:,1);
        %------------------------------------------------
        user_info.ids=obj.data.varobs_id; % location of the observables in the state vector
        locs=locate_variables(targets,obj.observables.name); % locations of the chosen variables
        user_info.ids=user_info.ids(locs); % locations of the chosen variables in the state vector
        y=data(locs,:);
        [ny,smpl]=size(y);
        if isempty(prior_sample)
            prior_sample=smpl;
        end
        user_info.T=prior_sample;
        
        zero_frequency_spectral_density();
        
        function zero_frequency_spectral_density()
            % demean the data and address problem of missing obs
            %----------------------------------------------------
            for ii=1:ny
                if ~obj.options.data_demean
                    di=y(ii,:);
                    mi=mean(di(~isnan(di)));
                    y(ii,:)=di-mi;
                end
                %  Address the problem of missing observations somehow. For
                %  the moment we just set missings to the mean i.e. 0
                %  (demeaned variables). An alternative could be to
                %  interpolate
                %  --------------------------------------------------------
                missing=isnan(y(ii,:));
                y(ii,missing)=0;
            end
            dyy=y.*y;
            % variances to match
            %--------------------
            user_info.Fhat=sum(dyy,2)/smpl;
            hF=hfunc(user_info.Fhat);
            % zero-frequency spectral density
            %--------------------------------
            user_info.Shat=autocovariance(0);
            for l=1:user_info.ll
                Cl=autocovariance(l);
                user_info.Shat=user_info.Shat+(1-l/(user_info.ll+1))^user_info.theta*(Cl+Cl');
            end
            function Cj=autocovariance(jj)
                Cj=1/(smpl-jj)*(hF(:,1:end-jj)*hF(:,jj+1:end)');
            end
            function h=hfunc(F)
                h=bsxfun(@minus,dyy,F);
            end
        end
    end
end
