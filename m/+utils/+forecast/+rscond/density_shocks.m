function [shocks,retcode]=density_shocks(R,mu_lb_ub,nshocks,OMG,options)
% Both gibbs and ghk algorithms seem to work fine regardless of whether
% there is shock simul_shock_uncertainty or not. But then again, only when the default
% algorithm, which does not trim the R (shock impact) matrix is employed.

[~,defaults]=utils.forecast.rscond.tmvnrnd();
defaults=[defaults
    {
    'simul_shock_uncertainty',false,@(x)islogical(x),...
    'simul_shock_uncertainty must be a logical'
    }
    ];

if nargin==0
    shocks=cell2struct(defaults(:,2),defaults(:,1),1);
    retcode=defaults;
    return
else
    narginchk(2,5)
end

if nargin<5
    options=[];
    if nargin<4
        OMG=[];
        if nargin<3
            nshocks=[];
        end
    end
end
if isempty(options)
    options=struct();
end

options=parse_arguments(defaults,options);
simul_shock_uncertainty=options.simul_shock_uncertainty;
options=rmfield(options,'simul_shock_uncertainty');

if isempty(OMG)
    OMG=R*R.';
end
mu=mu_lb_ub(:,1);
lb=mu_lb_ub(:,2);
ub=mu_lb_ub(:,3);

x=utils.forecast.rscond.tmvnrnd(mu,OMG,lb,ub,options);

[M1,M2,RM2i,~,~,~,retcode]=utils.forecast.rscond.null_and_column_spaces(R,options.debug);

if retcode
    shocks=[];
else
    M2_gam2=hard_conditions(x);
    
    shocks=M2_gam2;
    if simul_shock_uncertainty
        n1=size(M1,2);
        gam1=randn(n1,options.forecast_conditional_sampling_ndraws);
        shocks=M1*gam1+shocks;
    end
    
    if options.debug
        % all lb<=x<=ub
        %---------------
        fprintf('Checking that lb<=R*shocks: %0.0f\n\n',allpos(R*shocks,lb))
        
        fprintf('Checking that R*shocks<=ub : %0.0f\n\n',allpos(ub,R*shocks))
    end
    shocks=reshape_shocks(shocks);
end

    function shocks=reshape_shocks(A)
        if ~isempty(nshocks)
            [nr,np]=size(A);
            nc=nr/nshocks;
            shocks=reshape(A,[nshocks,nc,np]);
            if options.debug
                shocks2=zeros(nshocks,nc,np);
                for ipage=1:np
                    shocks2(:,:,ipage)=reshape(A(:,ipage),nshocks,nc);
                end
                shocks3=zeros(nshocks,nc,np);
                offset=0;
                for ipage=1:np
                    shocks3(offset+(1:nr))=A(:,ipage);
                    offset=offset+nr;
                end
                disp('checking reshape shocks')
                all(all(all(shocks-shocks2==0)))
                all(all(all(shocks-shocks3==0)))
            end
        end
    end

    function ehat=hard_conditions(yhats)
        ehat=M2*RM2i*yhats;
        % check that ehat is a solution to our problem
        %----------------------------------------------
        if options.debug
            fprintf('Hard-condition problem max resid = %0.8f\n\n',...
                maxdiff(R*ehat,yhats))
            
            % M2*gam2 and R'*inv(R*R')*x
            %------------------------------
            RRRi=utils.forecast.rscond.direct_inverse(R);
            ehat2=RRRi*yhats;
            fprintf('M2*gam2 vs R''*inv(R*R'')*x max diff = : %0.8f\n\n',...
                maxdiff(ehat,ehat2))
            
            % pinv(R)*x and R'*inv(R*R')*x
            %------------------------------
            ehat3=pinv(full(R))*x;
            fprintf('pinv(R)*x vs R''*inv(R*R'')*x max diff = %0.8f\n\n',...
                maxdiff(ehat3,ehat2))
            
            % pinv(R)*x and M2*gam2
            %-----------------------
            fprintf('pinv(R)*x vs M2*gam2 max diff = %0.8f\n\n',...
                maxdiff(ehat3,ehat))
            
            % R'*((R*R')\yhats) and R'*inv(R*R')*x
            %--------------------------------------
            ehat4=R'*((R*R')\yhats);
            fprintf('R''*((R*R'')\\x) vs R''*inv(R*R'')*x max diff = %0.8f\n\n',...
                maxdiff(ehat4,ehat2))
            
            % R'/(R*R')*yhats and R'*inv(R*R')*x
            %--------------------------------------
            ehat5=R'/(R*R')*yhats;
            fprintf('R''/(R*R'')*x vs R''*inv(R*R'')*x max diff = %0.8f\n\n',...
                maxdiff(ehat5,ehat2))
            
            % R'/(R*R')*yhats and pinv(R)
            %--------------------------------------
            fprintf('R''/(R*R'')*x vs pinv(R)*x max diff = %0.8f\n\n',...
                maxdiff(ehat5,ehat3))
        end
    end

    function flag=allpos(a,b)
        flag=all(all(abs(bsxfun(@minus,a,b))>=0));
    end

    function d=maxdiff(a,b)
        d=max(max(abs(bsxfun(@minus,a,b))));
        d=full(d);
    end
end