function obj=unclassified(obj,varargin)

if isempty(obj)
    
    mydefaults=the_defaults();
    
    if nargout
        
        obj=mydefaults;
        
    else
        
        clear obj
        
        disp_defaults(mydefaults);
        
    end
    
    return
    
end

end

function d=the_defaults()

num_fin=@(x)isnumeric(x) && isfinite(x);

num_fin_int=@(x)num_fin(x) && floor(x)==ceil(x) && x>=0;

rows=3;

cols=3;

d={
    'results_folder','',@(x)ischar(x),...
    'results_folder must be a char'
    
    'graphics',[rows,cols],...
    @(x) isequal(size(x),[1,2]) && all(x>0) && num_fin_int(x(1)) && ...
    num_fin_int(x(2)),...
    'graphics must be a 1 x 2 vector of integers'
    
    'debug(rs)',false,@(x)islogical(x),'debug must be a logical'
    };            

end
