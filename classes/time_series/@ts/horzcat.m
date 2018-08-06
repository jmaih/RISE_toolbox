function this=horzcat(varargin)
% Combines/concatenate databases into a combined database
%
% Note:
%    - It is assumed that the frequency and length of time series are the same between the two databases.
%


for ii=1:length(varargin)
    
    if ii==1
        
        this=varargin{ii};
        
        if ~((isa(this,'double') && isempty(this))||isa(this,'ts'))
            
            error('first argument must be an empty double or a ts')
            
        end
        
    else
        
        if ~isequal(class(varargin{ii}),'ts')
            
            error([mfilename,':: input ',int2str(ii),' must be from class ts'])
            
        end
        
        if ~isempty(varargin{ii})
            
            this=this&varargin{ii};
            
        end
        
    end
    
end

end
