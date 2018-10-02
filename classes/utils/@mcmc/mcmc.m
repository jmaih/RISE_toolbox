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
        
        function obj=mcmc(draws,pnames,drop,start_from,trimming,ilinres)
            % Constructor for mcmc object
            %
            % ::
            %
            %    mcmc_helper = mcmc(draws, pnames, drop, start_from, trimming)
            %
            % Args:
            %
            %    draws (struct): output from samplers
            %
            %    pnames (cellstr): cell of parameter names
            %
            %    drop (double): fraction of samples to drop (default: 0.5)
            %
            %    start_from (integer): discard first (start_from-1) samples (default: 1)
            %
            %    trimming (integer): only use every trimming value (default: 1)
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
            
            if n<6
                
                ilinres=[];
                
                if n<5
                    
                    trimming=[];
                    
                    if n<4
                        
                        start_from=[];
                        
                        if n<3
                            
                            drop=[];
                            
                            if n<2
                                
                                pnames=[];
                                
                                if n<1
                                    
                                    obj=mcmc.empty();
                                    
                                    return
                                    
                                end
                                
                            end
                            
                        end
                        
                    end
                    
                end
                
            end
            
            [obj.draws,obj.i_dropped,obj.start,results_summary]=mcmc.process_draws(draws,drop,start_from,trimming);
            
            obj.best=results_summary.best_of_the_best;
            
            if ~isempty(ilinres)
                % untransform the parameters
                for ii=1:numel(obj.draws)
                    
                    obj.draws(ii).x=ilinres(obj.draws(ii).x);
                    
                end
                
            end
            
            [obj.nchains,obj.npop]=size(obj.draws);
            
            if obj.npop
                
                obj.nparams=numel(obj.draws(1).x);
                
            else
                
                obj=mcmc.empty();
                
                return
                
            end
            
            if isempty(pnames)
                
                pnames=parser.create_state_list('p',obj.nparams);
                
            end
            
            obj.pnames=pnames;
            
            if obj.npop
                
                obj.psrf=gelman_rubin(obj,true);
                
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
        
    end
    
    methods(Access=private)
        
        varargout=load_draws(varargin)
        
        varargout=gelman_rubin(varargin)
        
    end
    
end