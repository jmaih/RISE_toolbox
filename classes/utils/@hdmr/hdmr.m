classdef hdmr
    % hdmr High dimensional model representation
    %
    % hdmr Methods:
    % ----------------
    %
    % estimate -
    % first_order_effect -
    % hdmr -   objective is either : f and theta or a function that will
    % metamodel -
    % plot_fit -
    % polynomial_evaluation -   later on, the function that normalizes could come in here so that the
    % polynomial_integration -   polynomial is of the form a0+a1*x+...+ar*x^r
    % polynomial_multiplication -   each polynomial is of the form a0+a1*x+...+ar*x^r
    %
    % hdmr Properties:
    % -------------------
    %
    % N -
    % Nobs -
    % n -
    % output_nbr -
    % theta -
    % theta_low -
    % theta_high -
    % g -
    % x -
    % expansion_order -
    % pol_max_order -
    % poly_coefs -
    % Indices -
    % coefficients -
    % aggregate -
    % f0 -
    % D -
    % sample_percentage -
    % optimal -
    % param_names -
    properties
        N
        Nobs
        n
        output_nbr
        theta
        theta_low
        theta_high
        g
        x
        expansion_order
        pol_max_order
        poly_coefs
        Indices
        coefficients
        aggregate
        f0
        D
        sample_percentage
        optimal=true
        param_names
    end
    methods
        %         function day = get.pol_max_order(obj),day = obj.pol_max_order;end
        %         %================
        %         function obj = set.pol_max_order(obj,neworder),obj.pol_max_order = neworder;end
        %         %================
        function obj=hdmr(objective,param_names,bounds,expansion_order,pol_max_order,...
                pol_optimal,sample_percentage)
            % objective is either : f and theta or a function that will
            % help calculate f and theta
            if nargin<1
                return
            end
            if nargin<7
                sample_percentage=[];
                if nargin<6
                    pol_optimal=[];
                    if nargin<5
                        pol_max_order=[];
                        if nargin<4
                            expansion_order=[];
                            if nargin<3
                                bounds=[];
                                if nargin<2
                                    param_names=[];
                                end
                            end
                        end
                    end
                end
            end
            lb=bounds(:,1);
            ub=bounds(:,2);
            if isempty(param_names)
                param_names=cellstr(num2str(1:numel(lb)));
            elseif ischar(param_names)
                param_names=cellstr(param_names);
            end
            if numel(param_names)~=numel(lb)
                error([mfilename,':: number of parameter names should match numel(lb)'])
            end
            if ~iscell(objective)
                error([mfilename,':: objective should be a cell {fhandle,bounds,N} or {f,theta}'])
            end
            % check the type of inputs
            sampled=numel(objective)==2;
            if ~sampled
                if numel(objective)~=3
                    error([mfilename,':: when the variables are not sampled, objective should be {fhandle,bounds,N} '])
                end
                fhandle=objective{1};
                lb=objective{2}(:,1);
                ub=objective{2}(:,2);
                nn=numel(lb);
                NN=objective{3};
                theta_=quasi_monte_carlo.sobol(nn,NN,lb,ub); % number of sub-intervals
                f_=fhandle(theta_);
                objective={f_,theta_};
            end
            %=================
            obj.param_names=param_names;
            obj.g=objective{1};
            obj.theta=objective{2};
            [obj.output_nbr,obj.N]=size(obj.g);
            [obj.n,c2]=size(obj.theta);
            if c2~=obj.N
                error([mfilename,':: number of parameter simulations inconsistent with number of columns of function values'])
            end
            %=================
            [obj.x,obj.theta_low,obj.theta_high]=theta_to_x(obj.theta,lb,ub);
            %
            if isempty(expansion_order),expansion_order=2;end
            if isempty(pol_max_order),pol_max_order=3;end
            if isempty(sample_percentage),sample_percentage=.75;end
            if isempty(pol_optimal),pol_optimal=true;end
            
            sample_percentage=abs(sample_percentage);
            while sample_percentage>1
                sample_percentage=0.1*sample_percentage;
            end
            obj.sample_percentage=sample_percentage;
            obj.optimal=pol_optimal;
            obj.Nobs=ceil(obj.sample_percentage*obj.N);
            obj.expansion_order=expansion_order;
            obj.pol_max_order=pol_max_order;
            obj.Indices=cell(1,obj.expansion_order);
            obj.Indices{1}=(1:obj.n)';
            for lev=2:obj.expansion_order
                previous=obj.Indices{lev-1};
                for k=1:size(previous,1)
                    last=previous(k,end);
                    batch=(last+1:obj.n)';
                    nb=numel(batch);
                    if nb
                        obj.Indices{lev}=[
                            obj.Indices{lev}
                            [previous(k*ones(nb,1),:),batch]
                            ];
                    end
                end
            end
        end
        %================
        function obj=estimate(obj,debug)
            if nargin==1
                debug=false;
            end
            % constant (Expectation) term
            obj.f0=mean(obj.g(:,1:obj.Nobs),2);
            % the total variance
            obj.D=1/obj.Nobs*sum(bsxfun(@minus,obj.g(:,1:obj.Nobs),obj.f0).^2,2);
            %
            obj.poly_coefs=orthonormal_polynomial(obj.pol_max_order,obj.x(:,1:obj.Nobs),obj.optimal,debug);
            
            % computation of coefficients
            obj.coefficients=struct('level',{},'indexes',{},'regime',{},'coef',{},'sensitivity',{});
            obj.aggregate=struct('level',{},'indexes',{},'sensitivity',{});
            iter=0;
            aggr_iter=0;
            POLS=ones(obj.Nobs,obj.pol_max_order,obj.n);
            for v=1:obj.n
                xp=obj.x(v,1:obj.Nobs)';
                POLS(:,2,v)=xp;
                for oo=2:obj.pol_max_order
                    POLS(:,oo+1,v)=POLS(:,oo,v).*xp;
                end
            end
            
            for lev=1:obj.expansion_order
                indexes=obj.Indices{lev};
                for jj=1:size(indexes,1)
                    ind=indexes(jj,:);
                    ncols=numel(ind);
                    % now get the grid with each dimension going from 1 to pol_max_order.
                    % This exactly where one could think about reducing the obj.poly_coefs
                    % order to, say l<pol_max_order as the number of factors (ncols)
                    % increases.
                    Regimes=utils.gridfuncs.mygrid(obj.pol_max_order*ones(ncols,1));
                    %         nregs=size(Regimes,1);
                    aggr=0;
                    for kk=1:size(Regimes,1)
                        f_phi=transpose(obj.g(:,1:obj.Nobs));
                        % so that the number of outputs is in columns and
                        % the number of simulations in rows
                        for cc=1:ncols
                            vc=ind(cc); % variable
                            pol_ord_c=Regimes(kk,cc); % polynomial order
                            % %                             koefs=obj.poly_coefs{pol_ord_c}(:,vc);
                            koefs=obj.poly_coefs{pol_ord_c};%(:,vc)
                            regressor_c=hdmr.polynomial_evaluation(koefs,POLS(:,:,vc));
                            f_phi=bsxfun(@times,f_phi,regressor_c);
                        end
                        iter=iter+1;
                        obj.coefficients(iter).level=lev;
                        obj.coefficients(iter).indexes=ind;
                        obj.coefficients(iter).regime=Regimes(kk,:);
                        obj.coefficients(iter).coef=transpose(mean(f_phi,1)); % mean over the rows
                        % the coefs are vectors, depending on the number of
                        % outputs
                        obj.coefficients(iter).sensitivity=obj.coefficients(iter).coef.^2./obj.D;
                        aggr=aggr+obj.coefficients(iter).sensitivity;
                    end
                    aggr_iter=aggr_iter+1;
                    obj.aggregate(aggr_iter).level=lev;
                    obj.aggregate(aggr_iter).indexes=ind;
                    obj.aggregate(aggr_iter).sensitivity=aggr;
                end
            end
        end
        %================
        function fit=metamodel(obj,theta)
            if nargin<2
                theta='insample';
            end
            if ischar(theta)
                if strcmp(theta,'insample')
                    theta=obj.theta(:,1:obj.Nobs);
                    ff=obj.g(:,1:obj.Nobs);
                elseif strcmp(theta,'outofsample')
                    theta=obj.theta(:,obj.Nobs+1:end);
                    ff=obj.g(:,obj.Nobs+1:end);
                else
                    error([mfilename,':: wrong flag'])
                end
            elseif iscell(theta)
                ff=theta{1};
                theta=theta{2};
            end
            Neff=size(theta,2);
            % normalize between 0 and 1
            xx=theta_to_x(theta,obj.theta_low,obj.theta_high);
            
            fit=struct();
            fit.f0=obj.f0;
            fit.h=obj.f0(:,ones(Neff,1));
            fit.ff=ff;
            iter=0;
            aggr_iter=0;
            fit.theta=theta;
            for lev=1:obj.expansion_order
                indexes=obj.Indices{lev};
                for jj=1:size(indexes,1)
                    ind=indexes(jj,:);
                    ncols=numel(ind);
                    % now get the grid with each dimension going from 1 to pol_max_order.
                    % This exactly where one could think about reducing the poly_coefs
                    % order to, say l<pol_max_order as the number of factors (ncols)
                    % increases.
                    Regimes=utils.gridfuncs.mygrid(obj.pol_max_order*ones(ncols,1));
                    %         nregs=size(Regimes,1);
                    f_ijklmn=0;
                    fijk_name='f';
                    for kk=1:size(Regimes,1)
                        phi_pqrst=1;
                        for cc=1:ncols
                            vc=ind(cc); % variable
                            if kk==1
                                fijk_name=[fijk_name,'_',int2str(vc)]; %#ok<AGROW>
                            end
                            pol_ord_c=Regimes(kk,cc); % polynomial order
                            % %                             koefs=obj.poly_coefs{pol_ord_c}(:,vc);
                            koefs=obj.poly_coefs{pol_ord_c};%(:,vc)
                            regressor_c=hdmr.polynomial_evaluation(koefs,xx(vc,:)');
                            phi_pqrst=phi_pqrst.*regressor_c;
                        end
                        iter=iter+1;
                        thecoefs=transpose(obj.coefficients(iter).coef);
                        f_ijklmn=f_ijklmn+bsxfun(@times,thecoefs,phi_pqrst);
                    end
                    f_ijklmn=transpose(f_ijklmn);
                    aggr_iter=aggr_iter+1;
                    fit.(fijk_name)=f_ijklmn;
                    fit.h=fit.h+f_ijklmn;
                end
            end
        end
        %================
        function plot_fit(obj,theta)
            fit=obj.metamodel(theta);
            for oo=1:obj.output_nbr
                figure('name',['model fit for output #',int2str(oo)])
                g_=fit.ff(oo,:);
                h=fit.h(oo,:);
                subplot(3,1,1)
                axis([min(g_(:)),max(g_(:)),min(h(:)),max(h(:))])
                hold on
                scatter(g_(:),h(:))
                xlabel('True model')
                ylabel('Metamodel')
                plot([min(g_(:)),max(g_(:))],[min(h(:)),max(h(:))],'color','r')
                hold off
                
                subplot(3,1,2)
                [fg,xg]=distributions.empirical_cdf(g_,min(g_),max(g_),100);
                [fh,xh]=distributions.empirical_cdf(h,min(g_),max(g_),100);
                plot(xg,fg,xh,fh)
                title('cdf')
                legend({'True model','Metamodel'})
                
                subplot(3,1,3)
                [fg,xg]=distributions.kernel_density(g_,min(g_),max(g_),'normal',100);
                [fh,xh]=distributions.kernel_density(h,min(g_),max(g_),'normal',100);
                plot(xg,fg,xh,fh)
                title('pdf')
                legend({'True model','Metamodel'})
                
            end
        end
        %================
        function first_order_effect(obj,theta,index,output_index)
            if nargin<4
                output_index=1;
            end
            if ischar(index)
                pname=index;
                index=locate_variables(pname,obj.param_names);
            elseif isnumeric(index)
                pname=obj.param_names{index};
            end
            if ~ismember(output_index,(1:obj.output_nbr))
                error([mfilename,':: output_index out of range. Must be in [1,',int2str(obj.output_nbr),']'])
            end
            fit=obj.metamodel(theta);
            %                 figure('name',['first-order effect of variable ',int2str(index),' for output #',int2str(output_index)])
            g_=fit.ff(output_index,:);
            scatter(fit.theta(index,:),g_)
            hold on
            [xi,tag]=sort(fit.theta(index,:));
            f_i=fit.(['f_',int2str(index)]);
            f_i=f_i(tag);
            plot(xi,obj.f0+f_i,'r','linewidth',1.5)
            xlabel(['Input: ',pname])
            legend({'f',['f0+f_',int2str(index)]},'interpreter','none')
            title(['S_',int2str(index),'=',num2str(obj.aggregate(index).sensitivity)])
        end
    end
    methods(Static)
        function f=polynomial_evaluation(p0,x)
            % later on, the function that normalizes could come in here so that the
            % normalization is done according to the hdmr_type of polynomial chosen.
            f=p0(1);
            power_done=min(size(x))>1; % <--- size(x,2)>1
            % the quad function may input a row vector
            if power_done
                xi=x(:,2);
            else
                xi=x;
            end
            for ii=2:numel(p0)
                f=f+p0(ii)*xi;
                if ii<numel(p0)
                    if power_done
                        xi=x(:,ii+1);
                    else
                        xi=xi.*x;
                    end
                end
            end
        end
        function coefs=polynomial_multiplication(p1,p2)
            % each polynomial is of the form a0+a1*x+...+ar*x^r
            
            % get the orders
            o1=numel(p1)-1;
            o2=numel(p2)-1;
            
            nterms=o1+o2+1;
            coefs=zeros(1,nterms);
            
            for i1=1:o1+1
                for i2=1:o2+1
                    order=(i1-1)+(i2-1);
                    coefs(order+1)=coefs(order+1)+p1(i1)*p2(i2);
                end
            end
        end
        function val=polynomial_integration(p,a,b)
            % polynomial is of the form a0+a1*x+...+ar*x^r
            % the integral is then a0*x+a1/2*x^2+...+ar/(r+1)*x^(r+1)
            order=numel(p)-1;
            pbar=p(:)./(1:order+1)';
            b_a=b-a;
            val=0;
            % exploit horner
            for oo=order+1:-1:1
                val=(val+pbar(oo))*b_a;
            end
        end
    end
end


