classdef rise_endo_priors
    properties(Constant)
        ll=2;
        theta=1;
    end
    properties
        Shat % zero-frequency spectral density
        Fhat % variances to match
        smpl % length of data
        T % size of the prior sample, representing the strength of the beliefs
        n % number of targets
%         y % demeaned and de-nan-ed data
        ids % locations of targets in state vector
        targets % target names
        state_size % total number of endogenous variables in the rise model
    end
    methods
        function obj=rise_endo_priors(this_rise,targets,T)
            if nargin<3
                T=[];
                if nargin<2
                    targets=[];
                    if nargin<1
                        return
                    end
                end
            end
            if isempty(this_rise)
                return
            end
            if isempty(targets)
                targets=this_rise.observables.name;
            elseif ischar(targets)
                targets=cellstr(targets);
            end
            obj.targets=targets;
            obj.ids=this_rise.varobs.id; % location of the observables in the state vector
            locs=locate_variables(targets,this_rise.observables.name); % locations of the chosen variables
            obj.ids=obj.ids(locs); % locations of the chosen variables in the state vector
            data=this_rise.data.y(:,this_rise.data.start:this_rise.data.finish);
            obj.state_size=this_rise.endogenous.number(2);
            
            obj=zero_frequency_spectral_density(obj,data,locs);
            if isempty(T)
                T=obj.smpl;
            end
            obj.T=T;
        end
        function ld=log_density(obj,P) % P is the theoretical variance of the model
            if ~isequal(size(P),[obj.state_size,obj.state_size])
                error([mfilename,':: wrong size of the covariance matrix. Expecting [',...
                    int2str(obj.state_size),',',int2str(obj.state_size),']'])
            end
            Ftheta=diag(P(obj.ids,obj.ids));
            FF=obj.Fhat-Ftheta;
            ld=.5*obj.n*log(obj.T/(2*pi))-.5*log(det(obj.Shat))-.5*obj.T*(FF'/obj.Shat)*FF;
        end
        function obj=zero_frequency_spectral_density(obj,data,locs)
            %=============demean the data==============
            y=data;
            for ii=1:size(data,1)
                di=data(ii,:);
                mi=mean(di(~isnan(di)));
                y(ii,:)=di-mi;
            end
            %==========================================
            y=y(locs,:);
            [obj.n,obj.smpl]=size(y);
            %  Now we really have to address the problem of missing observations somehow. The problem is
            %  especially prevalent when we mix data from multiple frequencies
            for ii=1:obj.n
                missing=isnan(y(ii,:));
                % just replace the missing by the uncoditional mean of the non-missing
                y(ii,missing)=mean(y(ii,~missing));
            end
            dyy=y.*y; % dyy=nan(n,smpl); for t=1:smpl, yt=y(:,t); dyy(:,t)=diag(yt*yt'); end            
            obj.Fhat=sum(dyy,2)/obj.smpl; % variance
            hF=hfunc(obj.Fhat);
            obj.Shat=autocovariance(0);
            for l=1:obj.ll
                Cl=autocovariance(l);
                obj.Shat=obj.Shat+(1-l/(obj.ll+1))^obj.theta*(Cl+Cl');
            end
            function Cj=autocovariance(jj)
                Cj=1/(obj.smpl-jj)*(hF(:,1:end-jj)*hF(:,jj+1:end)');
            end
            function h=hfunc(F)
                h=bsxfun(@minus,dyy,F);
            end
        end
    end
    
    methods(Static)        
    end
end