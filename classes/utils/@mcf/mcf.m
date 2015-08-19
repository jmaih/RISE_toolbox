classdef mcf < handle
    % mcf Monte Carlo Filtering
    %
    % mcf Methods:
    % ---------------
    %
    % addlistener - Add listener for event.
    % cdf -
    % cdf_plot -
    % correlation_patterns_plot -
    % delete - Delete a handle object.
    % eq -  == (EQ)   Test handle equality.
    % findobj - Find objects matching specified conditions.
    % findprop - Find property of MATLAB handle object.
    % ge -  >= (GE)   Greater than or equal relation for handles.
    % gt -  > (GT)   Greater than relation for handles.
    % isvalid - Test handle validity.
    % kolmogorov_smirnov_test -   tests the equality of two distributions using their CDFs
    % le -  <= (LE)   Less than or equal relation for handles.
    % lt -  < (LT)   Less than relation for handles.
    % mcf -   Example:
    % ne -  ~= (NE)   Not equal relation for handles.
    % notify - Notify listeners of event.
    % scatter -
    %
    % mcf Properties:
    % ------------------
    %
    % lb -   lower bounds for the parameters
    % ub -   upper bounds for the parameters
    % nsim -   number of simulations
    % procedure -   sampling procedure [{'uniform'}|'latin_hypercube'|'sobol'|'halton'|user-defined]
    % parameter_names -   names of the parameters
    % samples -
    % is_behaved -
    % nparam -
    % is_sampled -
    % check_behavior -
    % number_of_outputs -
    % user_outputs -
    % known_procedures -
    properties
        % lower bounds for the parameters
        lb
        % upper bounds for the parameters
        ub
        % number of simulations
        nsim=20000
        % sampling procedure [{'uniform'}|'latin_hypercube'|'sobol'|'halton'|user-defined]
        procedure
        % names of the parameters
        parameter_names
    end
    properties(SetAccess=protected)
        samples
        is_behaved
        nparam
        is_sampled = false
        check_behavior
        number_of_outputs
        user_outputs
    end
    properties(Constant)
        known_procedures={'uniform','latin_hypercube','sobol','halton'}
    end
    methods
        function obj=mcf(check_behavior,nsim_or_draws,lb,ub,names,procedure_)
            % Example:
            %   test=mcf({@(x)sum(x>0)==sum(x<=0),1},20000,-ones(10,1),ones(10,1))
            if nargin
                if nargin<6
                    procedure_=[];
                    if nargin<5
                        names=[];
                        if nargin<4
                            ub=[];
                            if nargin<3
                                lb=[];
                                if nargin<2
                                    error('insufficient number of arguments')
                                end
                            end
                        end
                    end
                end
                % objective function
                %--------------------
                nout=[];
                if iscell(check_behavior)
                    nout=check_behavior{2};
                    check_behavior=check_behavior{1};
                    assert(nout>0 && isfinite(nout) && ceil(nout)==floor(nout),'wrong specification of the number of output arguments')
                end
                if ~isa(check_behavior,'function_handle')
                    error('first argument must be a function handle')
                end
                if isempty(nout)
                    obj.number_of_outputs=nargout(check_behavior);
                    if obj.number_of_outputs<1
                        if obj.number_of_outputs<0
                            error('number of output arguments of the behavioral function should be finite and should map to a specific function')
                        else
                            error('number of output arguments of the behavioral function should be at least 1')
                        end
                    end
                else
                    obj.number_of_outputs=nout;
                end
                obj.check_behavior=check_behavior;
                % number of draws
                %----------------
                if isscalar(nsim_or_draws)
                    obj.nsim=nsim_or_draws;
                else
                    [obj.nparam,obj.nsim]=size(nsim_or_draws);
                    obj.samples=nsim_or_draws;
                    obj.is_sampled=true;
                end
                % lower and upper bounds on the parameters
                %-----------------------------------------
                if ~obj.is_sampled
                    if isempty(lb)||isempty(ub)
                        error('no samples provided and so lb and ub must be given')
                    end
                    [nrows,ncols]=size(lb);
                    if ~isequal(size(ub),[nrows,ncols])
                        error('upper and lower bounds must be of the same size')
                    end
                    if ~(all(isfinite(lb)) && all(isfinite(ub)))
                        error('lower and upper bounds must be finite')
                    end
                    if any(ub<lb)
                        error('ub cannot be lower than lb')
                    end
                    obj.lb=lb;
                    obj.ub=ub;
                    obj.nparam=nrows;
                end
                % parameter names
                %----------------
                if isempty(names)
                    names=strcat({'p_'},cellstr(num2str((1:obj.nparam)')));
                    names=cellfun(@(x)x(~isspace(x)),names,'uniformOutput',false);
                end
                if ischar(names),names=cellstr(names);end
                if ~iscellstr(names),error('names should be char or cellstr'),end
                if numel(names)~=obj.nparam
                    error('number of parameter names does not match the actual number of parameters')
                end
                obj.parameter_names=names(:).';
                if ~isempty(procedure_)
                    if isa(procedure_,'function_handle')
                        try
                            ntest=5;
                            test=procedure_(obj.lb,obj.ub,ntest);
                            assert(isequal(size(test),[obj.nparam,ntest]))
                        catch ME
                            error('A user-defined sampling procedure should accept 3 arguments which are lb, ub, nsim')
                        end
                    else
                        if ~ischar(procedure_)||any(strcmp(procedure_,{obj.known_procedures}))
                            disp(obj.known_procedures)
                            error('specification procedure expected to be one of the above')
                        end
                    end
                    obj.procedure=procedure_;
                else
                    obj.procedure = obj.known_procedures{1};
                end
            end
            estimate(obj);
        end
        function d=cdf(obj,pname)
            d=struct();
            p_id=strcmp(pname,obj.parameter_names);
            if ~any(p_id)
                error(['parameter "',pname,'" not found'])
            end
            d.x_good=transpose(obj.samples(p_id,obj.is_behaved));
            d.x_bad=transpose(obj.samples(p_id,~obj.is_behaved));
            d.lb=obj.lb(p_id);
            d.ub=obj.ub(p_id);
            [d.f_behave,d.x]=distributions.empirical_cdf(d.x_good,d.lb,d.ub);
            [d.f_non_behave]=distributions.empirical_cdf(d.x_bad,d.lb,d.ub);
        end
        function fig_handles=cdf_plot(obj,pnames,graph_nrows,graph_ncols,titel)
            if nargin<5
                titel=sprintf('%s :: Comparison of distributions',mfilename);
                if nargin<4
                    graph_ncols=3;
                    if nargin<3
                        graph_nrows=3;
                        if nargin<2
                            pnames=obj.parameter_names;
                        end
                    end
                end
            end
            fig_handles=utils.plot.multiple(@(xname)one_subplot(xname),...
                pnames,titel,graph_nrows,graph_ncols,...
                'FontSize',11,'FontWeight','normal');
            function [texname,legend_]=one_subplot(xname)
                % x-axis is the range (lb:ub)
                % y-axis is the cdf
                d=cdf(obj,xname);
                [pValue,Dn,xn,largest]=mcf.kolmogorov_smirnov_test(d);
                plot(d.x,[d.f_behave,d.f_non_behave],...
                    [xn,xn],[d.f_behave(largest),d.f_non_behave(largest)],'linewidth',2);
                axis([d.lb,d.ub,0,1])
                texname=[xname,' (Dn=',num2str(Dn),', P-value=',num2str(pValue),')'];
                legend_={'cdf behavior','cdf non-behavior','max dist.'};
            end
        end
        function hdl=correlation_patterns_plot(obj,names,type,pval_cutoff)
            if nargin<4
                pval_cutoff=[];
                if nargin<3
                    type=[];
                    if nargin<2
                        names=[];
                    end
                end
            end
            if isempty(pval_cutoff),pval_cutoff=.05;end
            if isempty(names),names=obj.parameter_names;end
            names_id=locate_variables(names,obj.parameter_names);
            if isempty(type),type='behave'; end
            if ~any(strcmp(type,{'behave','non-behave',''}))
                error('type must be one of the following: behave, non-behave, ''''')
            end
            switch type
                case 'behave'
                    data=obj.samples(names_id,obj.is_behaved);
                case 'non-behave'
                    data=obj.samples(names_id,~obj.is_behaved);
                otherwise
                    data=obj.samples(names_id,:);
            end
            npar=size(data,1);
            %---------------------------
            [corrmat,Pval]=corr(data.');
            hotties=Pval<=pval_cutoff;
            corrmat=tril(corrmat,-1);
            hotties(corrmat==0)=false;
%             [paired_names,~,pvalvec]=pairwise(names,corrmat,Pval);
            %----------------------------
            pax=nan(npar,npar);
            pax(hotties)=corrmat(hotties);
            
            titel=sprintf('%s :: Correlation patterns in the %s sample',mfilename,type);
            hdl=figure('name',titel);
            
            % do the plot
            %-------------
            imagesc(pax,[-1 1]);
            % remove label on y-axis
            %-----------------------
            tmp=gca();
            set(tmp,'YTickLabel','')
            
            % add the names of the parameters
            %--------------------------------
            for ip=1:npar
                text(ip,(0.5),names{ip},'HorizontalAlignment','left','rotation',90)%,'interpreter','none'
                text(0.5,ip,names{ip},'HorizontalAlignment','right')%,'interpreter','none'
            end
            
            % change the appearance
            %----------------------
            colorbar;
            ax=colormap;
            ax(1,:)=[0.9 0.9 0.9];
            colormap(ax);
            
            if npar>10
                set(tmp,'xtick',(5:5:npar))
                set(tmp,'ytick',(5:5:npar))
            end
            
            set(tmp,'xgrid','on')
            set(tmp,'ygrid','on')
        end
        function fig_handles=scatter(obj,names,type,pval_cutoff,graph_nrows,graph_ncols,titel)
            if nargin<7
                titel=[];
                if nargin<6
                    graph_ncols=3;
                    if nargin<5
                        graph_nrows=3;
                        if nargin<4
                            pval_cutoff=[];
                            if nargin<3
                                type=[];
                                if nargin<2
                                    names=[];
                                end
                            end
                        end
                    end
                end
            end
            if isempty(pval_cutoff),pval_cutoff=.05;end
            if isempty(names),names=obj.parameter_names;end
            names_id=locate_variables(names,obj.parameter_names);
            alternative=isempty(type);
            if isempty(titel)
                titel=sprintf('%s :: scatter plot of the data in %s the sample',mfilename,type);
            end
            data=obj.samples(names_id,:);
            if alternative
                select=obj.is_behaved;
            else
                switch type
                    case 'behave'
                        select=obj.is_behaved;
                    case 'non-behave'
                        select=~obj.is_behaved;
                    otherwise
                        error('expecting type to be "behave" or "non-behave"')
                end
            end
                names=names(:);
            % get all the correlations in the select sample and detect the
            % significant ones
            %-------------------------------------------------------------
            [corrmat,Pval]=corr(data(:,select).');
            [paired_names,~,pvalvec]=pairwise(names,corrmat,Pval);
            good=pvalvec<=pval_cutoff;
            if any(good)
                paired_names=paired_names(good);
                fig_handles=utils.plot.multiple(@(xname)do_one_scatter(xname),...
                    paired_names,titel,graph_nrows,graph_ncols,...
                    'FontSize',11,'FontWeight','normal');
            else
                warning(sprintf('no significant correlations found at %0.4f percent',100*pval_cutoff)) %#ok<SPWRN>
                fig_handles=[];
            end
            
            function [texname,legend_]=do_one_scatter(xname)
                arobase=find(xname=='@');
                pname1=xname(1:arobase-1);
                pname2=xname(arobase+1:end);
                p1=strcmp(pname1,names);
                p2=strcmp(pname2,names);
                do_the_scatter_plot(select,'b')
                legend_='';
                if alternative
                    hold on
                    do_the_scatter_plot(~select,'r')
                    legend_={'Behavior','Non-Behavior'};
                end
                axis tight
                xlabel(pname1,'interpreter','none')
                ylabel(pname2,'interpreter','none')
                texname=[num2str(corrmat(p1,p2)),'(',num2str(Pval(p1,p2)),')'];
                function do_the_scatter_plot(cols,couleur)
                    d1=data(p1,cols);
                    d2=data(p2,cols);
                    scatter(d1,d2,10,couleur,'filled','d')
                end
            end
        end
    end
    methods(Access=private)
        function estimate(obj)
            % draw all the parameters
            %-------------------------
            obj.get_draws();
            
            % initialize the waitbar
            %-----------------------
            x=struct('name','monte carlo filtering',...
                'message','initializing');
            utils.plot.waitbar('init',x)
            % outputs
            %----------
            obj.is_behaved=false(1,obj.nsim);
            nout=obj.number_of_outputs;
            output=cell(1,nout);
            obj.user_outputs=cell(obj.number_of_outputs,obj.nsim);
            for isim=1:obj.nsim
                draw=obj.samples(:,isim);
                [output{1:nout}]=obj.check_behavior(draw);
                obj.is_behaved(isim)=output{1};
                obj.user_outputs(:,isim)=output(:);
                utils.plot.waitbar('update',isim/obj.nsim)
            end
            utils.plot.waitbar('close')
            disp([num2str(100*sum(obj.is_behaved)/obj.nsim),' percent of the prior support is consistent with the behavior'])
        end
        function get_draws(obj)
            if ~obj.is_sampled
                if ischar(obj.procedure)
                    switch obj.procedure
                        case 'uniform'
                            LB=obj.lb(:,ones(obj.nsim,1));
                            UB=obj.ub(:,ones(obj.nsim,1));
                            obj.samples=LB+(UB-LB).*rand(obj.nparam,obj.nsim);
                        case {'latin_hypercube','sobol','halton'}
                            obj.samples=quasi_monte_carlo.(obj.procedure)(obj.lb,obj.ub,obj.nsim);
                        otherwise
                            error('unknown sampling procedure')
                    end
                else
                    obj.samples=obj.procedure(obj.lb,obj.ub,obj.nsim);
                end
            end
        end
    end
    methods(Static)
        function [pValue,Dn,xn,largest]=kolmogorov_smirnov_test(d)
            % tests the equality of two distributions using their CDFs
            ff=abs(d.f_behave-d.f_non_behave);
            largest=ff==max(ff);
            largest=find(largest);
            if numel(largest)>1
                disp([parname,' has more than one location with largest distance'])
                largest=largest(1);
            end
            Dn=ff(largest);
            xn=d.x(largest);
            %----------------
            n1     =  length(d.x_good);
            n2     =  length(d.x_bad);
            n      =  n1 * n2 /(n1 + n2);
            lambda =  max((sqrt(n) + 0.12 + 0.11/sqrt(n)) * Dn , 0);
            %
            %  Use the asymptotic Q-function to approximate the 2-sided P-value.
            %
            jjj       =  (1:101)';
            pValue  =  2 * sum((-1).^(jjj-1).*exp(-2*lambda*lambda*jjj.^2));
            pValue  =  min(max(pValue, 0), 1);
            %----------------
        end
    end
end

function varargout=pairwise(varargin)
nin=nargin;
nout=nargout;
if nout>nin
    error('# output arguments cannot exceed # input arguments')
end
pnames=varargin{1}(:);
if ~iscellstr(pnames)
    error('first input argument must be a cellstr')
end
varargout=cell(1,nout);
n=numel(pnames);
m=sum(1:n-1);
varargout(2:end)=repmat({zeros(m,1)},1,nout-1);
cnames=cell(m,1);
offset=0;
start=n-1;
for ii=1:n-1
    cnames(offset+(1:start))=strcat(pnames{ii},'@',pnames(ii+1:end));
    for iout=2:nout
        varargout{iout}(offset+(1:start))=varargin{iout}(ii+1:end,ii);
    end
    % hold ground
    %------------
    offset=offset+start;
    % next round
    %-----------
    start=start-1;
end
varargout{1}=cellfun(@(x)x(~isspace(x)),cnames,'uniformOutput',false);
end