function [x1,S,options]=tmvnrnd(m,SIG,lb,ub,options)

num_fin_int=@(x)isnumeric(x)&&isfinite(x)&&isreal(x)&&isscalar(x) && x>=0;

Stmp=@(x)isempty(x)||isstruct(x)||(ismatrix(x) && size(x,1)==size(x,2));

algorithms={'rejection','gibbs','ghk'};

defaults={
    'forecast_conditional_sampling_algo','rejection',@(x)ismember(x,algorithms),...
    ['algorithm must be one of ',cell2mat(strcat(algorithms,'|'))]
    
    'forecast_conditional_sampling_ndraws',1,@(x)num_fin_int(x) && x>0,...
	'forecast_conditional_sampling_ndraws must be a positive scalar'
    
    'forecast_conditional_sampling_burnin',0,@(x)num_fin_int(x),...
	'forecast_conditional_sampling_burnin must be a positive scalar'
    
    'forecast_conditional_sampling_thinning',1,@(x)num_fin_int(x) && x>0,...
	'forecast_conditional_sampling_thinning must be a positive scalar'
    
    'forecast_conditional_sampling_max_rejection_attempts',500,@(x)num_fin_int(x),...
	'forecast_conditional_sampling_max_rejection_attempts must be a positive scalar'
    
    'S',[],@(x)Stmp(x),'temporary variable S must be empty|struct|matrix'
    
    'debug',false,@(x)islogical(x),'debug must be logical'
    };

if nargin<5
    options=[];
    if nargin==0
        x1=cell2struct(defaults(:,2),defaults(:,1),1);
        S=defaults;
        return
    end
end
if isempty(options)
    options=struct();
end

options=parse_arguments(defaults,options);

forecast_conditional_sampling_max_rejection_attempts=options.forecast_conditional_sampling_max_rejection_attempts;

forecast_conditional_sampling_algo=options.forecast_conditional_sampling_algo;

forecast_conditional_sampling_ndraws=options.forecast_conditional_sampling_ndraws;

debug=options.debug;

forecast_conditional_sampling_thinning=options.forecast_conditional_sampling_thinning;

forecast_conditional_sampling_burnin=options.forecast_conditional_sampling_burnin;
if strcmp(forecast_conditional_sampling_algo,'rejection')
    forecast_conditional_sampling_burnin=0;
    forecast_conditional_sampling_thinning=1;
end

nrows=size(lb,1);

% use rejection if possible
%---------------------------
if all(lb==-inf) && all(ub==inf)
    if ~strcmp(forecast_conditional_sampling_algo,'rejection')
        warning('No constraints, switching to the more efficient "rejection" algorithm')
    end
    forecast_conditional_sampling_algo='rejection';
end
ishard=(ub-lb)<1e-8;
is_soft=~ishard;

x1=zeros(nrows,forecast_conditional_sampling_ndraws);
% hard conditions
%----------------
x1(ishard,:)=do_hard_conditions(lb(ishard));

nsoft=sum(is_soft);

% quick exit
%-----------
if nsoft==0
    return
end

S=options.S;
% soft conditions
%----------------
x1(is_soft,:)=do_soft_conditions(m(is_soft),lb(is_soft),ub(is_soft),...
    SIG(is_soft,is_soft));

    function xs=do_soft_conditions(ms,lbs,ubs,SIGs)
        if isempty(S)
            S=chol(SIGs,'lower');
        end
        x=ms;
        xs=zeros(nsoft,forecast_conditional_sampling_ndraws);
        start=1-forecast_conditional_sampling_burnin;
        for idraw=start:forecast_conditional_sampling_ndraws*forecast_conditional_sampling_thinning
            switch forecast_conditional_sampling_algo
                case 'rejection'
                    do_rejection()
                case 'gibbs'
                    do_gibbs()
                case 'ghk'
                    do_ghk()
                otherwise
                    error(['unknown algorithm ',forecast_conditional_sampling_algo])
            end
            if idraw>0 && rem(idraw,forecast_conditional_sampling_thinning)==0
                xs(:,idraw/forecast_conditional_sampling_thinning)=x;
            end
        end
        
        function do_ghk()
            eta=nan(nsoft,1);
            for ii=1:nsoft
                mi=ms(ii)+S(ii,1:ii-1)*eta(1:ii-1);
                ai=normalize(lbs(ii));
                bi=normalize(ubs(ii));
                eta(ii)=utils.forecast.rscond.tnrnd(ai,bi);
            end
            x=ms+S*eta;
            function an=normalize(a)
                an=(a-mi)/S(ii,ii);
            end
        end
        
        function do_rejection
            iter=0;
            success=false;
            while ~success && iter<forecast_conditional_sampling_max_rejection_attempts
                iter=iter+1;
                x=ms+S*randn(nsoft,1);
                success=all(x>=lbs & x<=ubs);
                if success
                    if debug
                        fprintf('Rejection: success after %0.0f attempts\n',iter)
                    end
                end
            end
            if debug && ~success
                fprintf('Rejection: no success after %0.0f attempts\n',forecast_conditional_sampling_max_rejection_attempts)
            end
        end
        
        function do_gibbs()
            % Jayesh Hukumchand Kotecha and Petar M. Djuric (1999): "Gibbs
            % sampling approach for generation of truncated multivariate
            % Gaussian random variables" Proceedings of the 1999 IEEE
            % International Conference on Acoustics, Speech, and Signal
            % Processing, pp:1757-1760
            for ii=1:nsoft
                mi=conditional_mean();
                x(ii)=utils.forecast.rscond.tnrnd(lbs(ii),ubs(ii),mi,S.si{ii});
                if options.debug && (x(ii)<lbs(ii)||x(ii)>ubs(ii))
                    keyboard
                end
            end
            function mi=conditional_mean()
                jj=[1:ii-1,ii+1:nrows];
                Cij=SIG(ii,jj);
                Cji=SIG(jj,ii);
                Cii=SIG(ii,ii);
                if idraw==start
                    if ii==1
                        S=struct('si',{cell(1,nrows)},...
                            'Cij_iCjj',{cell(1,nrows)});
                    end
                    iCjj=SIG(jj,jj)\eye(nrows-1);
                    S.Cij_iCjj{ii}=Cij*iCjj;
                    Vii=Cii-S.Cij_iCjj{ii}*Cji;
                    if Vii<0
                        warning(['negative conditional variance ',...
                            num2str(Vii),' for item ',int2str(ii),...
                            ' set to zero'])
                        S.si{ii}=0;
                    else
                        S.si{ii}=sqrt(Vii);
                    end
                end
                mi=ms(ii)+S.Cij_iCjj{ii}*(x(jj)-ms(jj));
            end
        end
    end
    function xh=do_hard_conditions(lbh)
        xh=lbh(:,ones(1,forecast_conditional_sampling_ndraws));
    end
end
