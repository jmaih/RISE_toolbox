classdef mcf < handle
% MCF Create a Monte Carlo Filtering
%
% mcf Methods:
% ---------------
%
% cdf - cumulative distribution function
% cdf_plot - plot of cdf
% correlation_patterns_plot - plot of correlation patterns
% kolmogorov_smirnov_test - test of equality of distributions
% mcf - creates an mcf object
% scatter - scatter plot of the data
%
% mcf Properties:
% ------------------
%
% lb -   lower bounds for the parameters
% ub -   upper bounds for the parameters
% nsim -   number of simulations
% procedure -   sampling procedure [{'uniform'}|'latin_hypercube'|'sobol'|'halton'|user-defined]
% parameter_names -   names of the parameters
% samples - parameter draws
% is_behaved - boolean flag for behaved parameter vectors
% nparam - number of parameters
% is_sampled - true if draws are available
% check_behavior - checks whether the vectors should be true or false
% number_of_outputs - number of outputs to check_behavior
% user_outputs - sampled additional outputs
% known_procedures - known sampling procedures
properties

graph_nrows=3

graph_ncols=3

end

properties(SetAccess=protected)
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
% sampled draws
samples
% Boolean vector describing the draws
is_behaved
% number of parameters
nparam
% number of output of the behavior function
number_of_outputs
% user requested additional outputs
user_outputs
end

properties(SetAccess=protected,Hidden)

is_sampled = false

check_behavior

end

properties(Constant)
% Default sampling procedures: 'uniform', 'latin_hypercube',
% 'sobol', 'halton'
known_procedures={'uniform','latin_hypercube','sobol','halton'}

end

methods

function obj=mcf(check_behavior,nsim_or_draws,lb,ub,names,procedure_)
% MCF -- Monte Carlo Filtering
%
% ::
%
%
%   obj=MCF(check_behavior,nsim_or_draws)
%
%   obj=MCF(check_behavior,nsim_or_draws,lb)
%
%   obj=MCF(check_behavior,nsim_or_draws,lb,ub)
%
%   obj=MCF(check_behavior,nsim_or_draws,lb,ub,names)
%
%   obj=MCF(check_behavior,nsim_or_draws,lb,ub,names,procedure_)
%
% Args:
%              %
%              % - **check_behavior** [function_handle|cell|vector]: (i) when
%              % it is a function handle, the function takes as input a vector
%              % of parameters and returns in its first output a boolean that
%              % is true if the paramter satisfies the behavior and false
%              % otherwise. The procedure is going to check the number of
%              % output arguments that the function returns and collect all
%              % those extra output arguments while sampling. (ii) when it is
%              % a cell, it should be a two-element cell such that the first
%              % element is the function handle and the second is the number
%              % of outputs desired by the user. (iii) if it is a vector, all
%              % elements are either 0 or 1 or boolean, describing whether
%              % each parameter vector checks the behavior or not.
%              %
%              % - **nsim_or_draws** [integer|matrix]: When it is an integer,
%              % nsim_or_draws is the number of draws to sample. When it is a
%              % matrix, it is the draws. No further sampling will be
%              % performed.
%              %
%              % - **lb** [empty|vector]: lower bound of the search space
%              %
%              % - **ub** [empty|vector]: upper bound of the search space
%              %
%              % - **names** [empty|char|cellstr]: names of the parameters. If
%              % empty, the names are created as p_i, where "i" is the order
%              % of the parameter in the list.
%              %
%              % - **procedure_** [{'uniform'}|'latin_hypercube'|'sobol'|...
%              % 'halton'|user-defined]: when it is user-defined, it should be
%              % a function handle, and should take 3 inputs(lb, ub, nsim) and
%              % return a parameter draw.
%              %
% Returns:
%    :
%              %
%              % - **obj** [mcf object]: containing, among other things, the
%              % samples, a flag determining their behavior.
%              %
% Note:
%              %
% Example:
%
% See also:

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

check_behave_error='check_behavior must be a function handle, a vector of zeros and ones, or a vector of logicals';

if iscell(check_behavior)

nout=check_behavior{2};

check_behavior=check_behavior{1};

assert(nout>0 && isfinite(nout) && ceil(nout)==floor(nout),'wrong specification of the number of output arguments')

end

if isnumeric(check_behavior)

if all(check_behavior==1|check_behavior==0)

check_behavior=logical(check_behavior);

else

error(check_behave_error)

end

end

is_already_checked=islogical(check_behavior);

if is_already_checked

nout=0;

else

if ~isa(check_behavior,'function_handle')

error(check_behave_error)

end

end

if isempty(nout)

obj.number_of_outputs=nargout(check_behavior);

if obj.number_of_outputs<1

nout=0;

for ii=1:4

myout=cell(1,ii);

try

[myout{1:ii}]=check_behavior(0.5*(lb+ub)); %#ok<NASGU>

nout=nout+1;

catch

break

end

end

if nout==0

error('not enough output arguments for the behavioral function')

end

obj.number_of_outputs=nout;

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

obj.nsim=size(nsim_or_draws,2);

end

if is_already_checked && obj.is_sampled

if obj.nsim~=numel(check_behavior)

error('number of elements in check_behavior does not match the sample size')

end

obj.is_behaved=check_behavior(:).';

obj.check_behavior=[];

else

end

% lower and upper bounds on the parameters
%-----------------------------------------
if obj.is_sampled

if isempty(lb)||isempty(ub)

lb=min(obj.samples,[],2);

ub=max(obj.samples,[],2);

end

nrows=numel(lb);

else

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

end

obj.nparam=nrows;

obj.lb=lb;

obj.ub=ub;
% create parameter names if empty!!!
%-----------------------------------
names=utils.char.create_names(names,'p',obj.nparam);

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
% cdf -- Cumulative distribution function
%
% ::
%
%
%   d=cdf(obj,pname)
%
% Args:
%              %
%              % - **obj** [mcf object]:
%              %
%              % - **pname** [char]: name of the parameter which one wants
%              % the cdf for.
%              %
% Returns:
%    :
%              %
%              % - **d** [struct]: with fields lb, ub, x, f_behave,
%              % f_non_behave, x_good, x_bad
%              %
% Note:
%              %
% Example:
            %
            % See also:
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
        
        function hdl=cdf_plot(obj,pnames,titel)
            % cdf_plot -- plot of the cdfs
            %
            % Syntax
            % -------
            % ::
            %
            %   hdl=cdf_plot(obj)
            %
            %   hdl=cdf_plot(obj,pnames)
            %
            %   hdl=cdf_plot(obj,pnames,titel)
            %
            % Inputs
            % -------
            %
            % - **obj** [mcf object]: 
            %
            % - **pnames** [empty|char|cellstr]: names of the parameters of
            % interest
            %
            % - **titel** [empty|char]: title of the figures 
            %
            % Outputs
            % --------
            %
            % - **hdl** [vector]: handles to the plotted figures
            % 
            % More About
            % ------------
            %
            % Examples
            % ---------
            %
            % See also:
            if nargin<3
                
                titel=sprintf('%s :: Comparison of distributions',mfilename);
                
                if nargin<2
                
                    pnames=obj.parameter_names;
                
                end
                
            end
            
            if ischar(pnames),pnames=cellstr(pnames); end
            
            if numel(pnames)>1
                
            hdl=utils.plot.multiple(@(xname)one_subplot(xname),...
                pnames,titel,obj.graph_nrows,obj.graph_ncols,...
                'FontSize',11,'FontWeight','normal');
            
            else
                
                if iscellstr(pnames)
                    
                    pnames=pnames{1};
                    
                end
                
                [~,~,hdl]=one_subplot(pnames);
                
                title(pnames)
                
            end
            
            function [texname,legend_,hdl]=one_subplot(xname)
                % x-axis is the range (lb:ub)
                % y-axis is the cdf
                d=cdf(obj,xname);
            
                [pValue,Dn,xn,largest]=mcf.kolmogorov_smirnov_test(d);
                
                hdl=plot(d.x,[d.f_behave,d.f_non_behave],...
                    [xn,xn],[d.f_behave(largest),d.f_non_behave(largest)],'linewidth',2);
                
                axis([d.lb,d.ub,0,1])
                
                texname={xname;[' (Dn=',num2str(Dn),', P-value=',num2str(pValue),')']};
                
                legend_={'cdf behavior','cdf non-behavior','max dist.'};
            
            end
            
        end
        
        function hdl=correlation_patterns_plot(obj,names,type,pval_cutoff)
            % correlation_patterns_plot -- plot of correlation patterns
            %
            % Syntax
            % -------
            % ::
            %
            %   hdl=correlation_patterns_plot(obj)
            %
            %   hdl=correlation_patterns_plot(obj,names)
            %
            %   hdl=correlation_patterns_plot(obj,names,type)
            %
            %   hdl=correlation_patterns_plot(obj,names,type,pval_cutoff)
            %
            % Inputs
            % -------
            %
            % - **obj** [mcf object]: 
            %
            % - **names** [empty|char|cellstr]: names of the parameters of
            % interest
            %
            % - **type** [empty|'behave'|'non-behave']: if empty, all the
            % sample is considered. If 'behave' only the behavior sample is
            % considered. If 'non-behave', only the non-behavior sample is
            % considered.
            %
            % - **pval_cutoff** [numeric|{0.05}]: cutoff for significant
            % correlations
            %
            % Outputs
            % --------
            %
            % - **hdl** [vector]: handles to the plotted figures
            % 
            % More About
            % ------------
            %
            % Examples
            % ---------
            %
            % See also:
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
            
            if isempty(type)
                
                type='';
                
            else
                
                if ~any(strcmp(type,{'behave','non-behave'}))
                    
                    error('type must be one of the following: behave, non-behave, ''''')
                    
                end
                
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
        
        function hdl=scatter(obj,names,type,pval_cutoff,titel)
            % scatter -- scatter plot of the data
            %
            % Syntax
            % -------
            % ::
            %
            %   hdl=scatter(obj)
            %
            %   hdl=scatter(obj,names)
            %
            %   hdl=scatter(obj,names,type)
            %
            %   hdl=scatter(obj,names,type,pval_cutoff)
            %
            %   hdl=scatter(obj,names,type,pval_cutoff,titel)
            %
            % Inputs
            % -------
            %
            % - **obj** [mcf object]: 
            %
            % - **names** [empty|char|cellstr]: names of the parameters of
            % interest
            %
            % - **type** [empty|'behave'|'non-behave']: if empty, all the
            % sample is considered. If 'behave' only the behavior sample is
            % considered. If 'non-behave', only the non-behavior sample is
            % considered.
            %
            % - **pval_cutoff** [numeric|{0.05}]: cutoff for significant
            % correlations
            %
            % - **titel** [empty|char]: title of the figures 
            %
            % Outputs
            % --------
            %
            % - **hdl** [vector]: handles to the plotted figures
            % 
            % More About
            % ------------
            %
            % Examples
            % ---------
            %
            % See also:
            if nargin<5
                
                titel=[];
                
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
                
                hdl=utils.plot.multiple(@(xname)do_one_scatter(xname),...
                    paired_names,titel,obj.graph_nrows,obj.graph_ncols,...
                    'FontSize',11,'FontWeight','normal');
            
            else
                
                warning(sprintf('no significant correlations found at %0.4f percent',100*pval_cutoff)) %#ok<SPWRN>
                
                hdl=[];
            
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
            
            if isempty(obj.is_behaved)
                % outputs
                %----------
                nout=obj.number_of_outputs;
                
                user_outputs_=cell(obj.number_of_outputs,obj.nsim);
                % embarassingly parallel
                %-----------------------
                NumWorkers=utils.parallel.get_number_of_workers();
                % initialize the waitbar
                %-----------------------
                if NumWorkers==0
                
                    x=struct('name','monte carlo filtering','message',...
                        'initializing');
                    
                    utils.plot.waitbar('init',x)
                
                end
                
                samples_=obj.samples;
                
                obj.samples=[];
                
                check_behavior_=obj.check_behavior;
                
                nsim_=obj.nsim;
                
                parfor (isim=1:obj.nsim,NumWorkers) % for isim=1:obj.nsim 
                
                    draw=samples_(:,isim);
                    
                    output=cell(1,nout);
                    
                    [output{1:nout}]=check_behavior_(draw); %#ok<PFBNS>
                    
                    user_outputs_(:,isim)=output(:);
                    
                    if NumWorkers==0
                    
                        utils.plot.waitbar('update',isim/nsim_)
                    
                    end
                    
                end
                
                obj.samples=samples_;
                
                obj.is_behaved=cell2mat(user_outputs_(1,:));
                
                obj.user_outputs=user_outputs_;
                
                if NumWorkers==0
                
                    utils.plot.waitbar('close')
                
                end
                
            end
            
            disp([num2str(100*sum(obj.is_behaved)/obj.nsim),' percent of the support is consistent with the behavior'])
        
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
%                 disp([parname,' has more than one location with largest distance'])
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