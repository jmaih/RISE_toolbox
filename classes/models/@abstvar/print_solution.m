function print_solution(obj,varlist)

if nargin<2
    
    varlist=[];
    
end

string='';

nobj=numel(obj);

for iobj=1:nobj
    
    if nobj>1
        
        string=int2str(iobj);
        
    end
    
    fprintf(1,'\n%s\n',['MODEL ',string,' SOLUTION']);
    
    print_low_level(obj(iobj),mfilename,varlist)
    
end

end
