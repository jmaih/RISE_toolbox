function this=setOptions(this,varargin)

n=length(varargin);

if n==0
    
    return
    
end

priorList={'minnesota','jeffrey',...
    'nw','normal-wishart','inw','sz'};

if rem(n,2)
    
    error('number of arguments must be pair')
    
end

for ii=1:2:n
    
    prop=varargin{ii};
    
    check_property(prop,varargin{ii+1});
    
end

    function v=check_property(p,v)
        
        switch p
            
            case 'debug'
                
                if ~(isscalar(v) && islogical(v))
                    
                    error('"debug" Must be a logical')
                    
                end
                
                this.(prop)=v;
                
            case 'data'
                
                if isa(v,'ts')
                    
                    v=pages2struct(v);
                    
                end
                
                check_data(v)
                
                this.(prop)=v;
                
%             case 'prior'
%                 
%                 check_prior(v)
%                 
%                 this.(prop)=v;
%                 
%                 
            case 'linear_restrictions'
                
                if ischar(v),v=cellstr(v); end
                
                if ~iscellstr(v)
                    
                    error('linear restrictions should be char or cellstr')
                    
                end
                
                inter=intersect(this.linear_restrictions_prime_time,v);
                
                if ~isempty(inter)
                    
                    disp(inter)
                    
                    error('The restrictions above cannot be set by the user')
                    
                end
                
                this.(prop)=[v(:).',this.linear_restrictions_prime_time];
                
            otherwise
                
                error(['unrecognized property ',p])
                
        end
        
    end

    function check_prior(v)
        
        if isempty(v)
            
            return
            
        end
        
        if ~isstruct(v)
            
            error('v must be a structure')
            
        end
        
        fd=fieldnames(v);
        
        if ~all(ismember(fd,fieldnames(vartools.prior_hyperparams)))
            
            error('unrecognized field(s) in prior')
            
        end
        
        if ~ismember(lower(v.type),priorList)
            
            error(['unrecognized type of prior type ',v.type])
            
        end
        
        if any(strcmpi(v.type,{'nw','normal-wishart'}))
            
            if ~isfield(v,'normal_wishart_eta')||isempty(v.normal_wishart_eta)
                
                error('wrong format for normal_wishart_eta')
                
            end
            
        elseif any(strcmpi(v.type,{'indep-normal-wishart','inw'}))
            
            if ~isfield(v,'independent_normal_wishart_eta')||isempty(v.independent_normal_wishart_eta)
                
                error('wrong format for independent_normal_wishart_eta')
                
            end
            
        end
        
    end

    function check_data(d)
        
        if ~isstruct(d)
            
            error('data must come in the form of a structure')
            
        end
        
        fd=fieldnames(d);
        
        if isempty(this.members)
            
            check_endogenous(fd)
            
        else
            
            if ~all(ismember(this.members,fd))
                
                error('all elements in members should be represented in the data')
                
            end
            
            for ig=1:numel(fd)
                
                check_endogenous(fieldnames(d.(fd{ig})),d.(fd{ig}))
                
            end
            
        end
        
        function check_endogenous(v,root)
            
            if nargin<2
                
                root=[];
                
            end
            
            if ~all(ismember(this.endogenous,v))
                
                error('all endogenous variables should have data')
                
            end
            
            if ~isempty(root)
                
                for iv=1:numel(this.endogenous)
                    
                    if ~isa(root.(this.endogenous{iv}),'ts')
                        
                        error('data must be time series')
                        
                    end
                    
                end
                
            end
            
        end
        
    end

end
