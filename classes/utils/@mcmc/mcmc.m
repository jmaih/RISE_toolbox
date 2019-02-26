classdef mcmc < handle
    % mcmc object makes diagnostic plots from mcmc draws
    %
    
    properties
        pnames
        nchains
        npop
        nparams
        draws
        i_dropped
        start
        psrf
        best
    end
    
    methods
        
        function obj=mcmc(draws,pnames,subset,ilinres)
            % Constructor for mcmc object
            %
            % ::
            %
            %    mcmc_helper = mcmc(draws)
            %    mcmc_helper = mcmc(draws, pnames)
            %    mcmc_helper = mcmc(draws, pnames, subset)
            %    mcmc_helper = mcmc(draws, pnames, subset, ilinres)
            %
            % Args:
            %
            %    draws (struct): output from samplers
            %
            %    pnames (cellstr): cell of parameter names
            %
            %    subset (cell array|{empty}): When not empty, subset is a
            %     1 x 2 cell array in which the first cell contains a
            %     vector selecting the columns to retain in each chain and
            %     the second column contains the chains retained. Any or
            %     both of those cell array containts can be empty. Whenever
            %     an entry is empty, all the information available is
            %     selected. E.g. subsetting with dropping and trimming
            %     mysubs={a:b:c,[1,3,5]}. In this example, the first
            %     element selected is the one in position "a" and
            %     thereafter every "b" element is selected until we reach
            %     element in position "c". At the same time, we select
            %     markov chains 1,3 and 5.
            %
            %    ilinres (function handle|{empty}): function handle that
            %       untransforms the parameters in the presence of linear
            %       restrictions.
            %
            % Returns:
            %    :
            %
            %    - **mcmc_helper** : mcmc object
            %
            % Note:
            %    - It is the responsibility of the user to provide the
            %      names of the parameters as this routine aims to be
            %      independent from any estimation procedure or class. If
            %      the priors are set in separate structure, their names
            %      can easily be obtained via ::
            %
            %         pnames = fieldnames(priors);
            %
            %    - Alternatively, if using a RISE object, parameter names
            %      can be obtained via::
            %
            %         pnames = model.estimation.priors.name;
            %
            %    - Note that burn-in option in samplers already discard
            %      values, so make sure that **start_from** parameters is the
            %      intended value.
            %
            
            n=nargin;
            
            if n<4
                
                ilinres=[];
                
                if n<3
                    
                    subset=[];
                    
                    if n<2
                        
                        pnames=[];
                        
                        if n<1
                            
                            obj=mcmc.empty();
                            
                            return
                            
                        end
                        
                    end
                    
                end
                
            end
            
            [obj.draws,~,results_summary]=mcmc.process_draws(draws,subset);
            
            obj.best=results_summary(1).best_of_the_best;
            
            for ii=2:numel(results_summary)
                
                if results_summary(ii).best_of_the_best.f<obj.best.f
                    
                    obj.best=results_summary(ii).best_of_the_best;
                    
                end
                
            end
            
            [~,obj.npop,obj.nchains]=size(obj.draws);
            
            if ~isempty(ilinres)
                % untransform the parameters
                
                for ipop=1:obj.npop
                    
                    for ic=1:obj.nchains
                        
                        dii=ilinres(obj.draws(:,ipop,ic));
                        
                        if ipop==1 && ic==1
                            
                            mydraws=dii(:,ones(1,obj.npop),ones(1,obj.nchains));
                            
                        end
                        
                        mydraws(:,ipop,ic)=dii;
                        
                    end
                    
                end
                
                obj.draws=mydraws; 
                
                clear mydraws
                
            end
            
            if obj.npop
                
                obj.nparams=size(obj.draws,1);
                
            else
                
                obj=mcmc.empty();
                
                return
                
            end
            
            if isempty(pnames)
                
                pnames=parser.create_state_list('p',obj.nparams);
                
            end
            
            obj.pnames=pnames;
            
            if obj.npop
                
                obj.psrf=gelman_rubin(obj);
                
            end
            
        end
        
        varargout=traceplot(varargin)
        
        varargout=densplot(varargin)
        
        varargout=meanplot(varargin)
        
        varargout=autocorrplot(varargin)
        
        varargout=psrf_plot(varargin)
        
        varargout=scatterplot(varargin)
        
        varargout=summary(varargin)
        
        % varargout=heidelberg_welch(varargin)
        
    end
    
    methods(Static)
        
        varargout=process_draws(varargin)
        
        varargout=reload_draws(varargin)
        
    end
    
    methods(Access=private)
        
        varargout=load_draws(varargin)
        
        varargout=gelman_rubin(varargin)
        
    end
    
end