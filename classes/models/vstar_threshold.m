classdef vstar_threshold < handle
    
    properties
        items=struct('transition_variable',{},'transition_description',{},...
            'type',{},'controlled',{},'threshold_priors',{});
    end
    
    properties(Constant)
        
        threshold_list={'exponential','logistic','logisticn','logistic2'}
        
    end
    
    methods
        
        function obj=vstar_threshold()
            
        end
        
        %         function disp(obj)
        %
        %             for it=1:numel(obj.items)
        %
        %                 disp(obj.items(it))
        %
        %             end
        %
        %         end
        
        function new(obj,tname,ttype,cpriors,tvars)
            
            if nargin<5
                    
                    tvars={};
                
                if nargin<4
                
                cpriors=[];
                    
                end
                
            end
            
            if isempty(tvars)
                
                tvars={};
                
            end
            
            if ischar(tvars)
                
                tvars=cellstr(tvars);
                
            end
            
            if ~isempty(tvars) && ~iscellstr(tvars)
                
                error('list must be char or cellstr')
                
            end
            
            nthresh=1;
            
            tdescript=tname;
            
            process_threshold()
            
            obj.items(end+1)=struct('transition_variable',tname,...
                'transition_description',tdescript,'type',...
                ttype,'controlled',{tvars},'threshold_priors',cpriors);
            
            function process_threshold()
                
                if iscell(ttype)
                    
                    if numel(ttype)~=2
                        
                        error('when cell, expecting two elements')
                        
                    end
                    
                    nthresh=ttype{2};
                    
                    ttype=ttype{1};
                    
                end
                
                if ~ischar(ttype)||~ismember(ttype,obj.threshold_list)
                    
                    disp(obj.threshold_list)
                    
                    error('Threshold type must be one of the above')
                    
                end
                
                switch ttype
                    
                    case {'exponential','logistic'}
                        
                        if nthresh~=1
                            
                            error([ttype,' allows for one threshold only'])
                            
                        end
                    
                    case 'logistic2'
                        
                        nthresh=2;
                        
                    case 'logisticn'
                        
                        if nthresh<1
                            
                            error([ttype,' allows for one or more thresholds'])
                            
                        end
                        
                end
                
                if isempty(cpriors)
                    
                    cpriors=nan(nthresh,2);
                    
                end
                
                if ~isequal(size(cpriors),[nthresh,2])
                    
                    error('size of transition thresholds inconsistent with transition function')
                    
                end
                
                cpriors=sortrows(cpriors,1);
                
                if iscell(tname)
                    
                    if numel(tname)~=2
                        
                        error('when cell, tname should have two components')
                        
                    end
                    
                    tdescript=tname{2};
                    
                    tname=tname{1};
                    
                end
                
            end
            
        end
        
    end
    
end