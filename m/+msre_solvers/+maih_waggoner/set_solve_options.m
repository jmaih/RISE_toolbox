function slvOpts=set_solve_options(varargin)
% INTERNAL FUNCTION
%

% [refine,checkStab,allSols,msvOnly,xplosRoots,debug]

defaults={'refine',false;
    'checkStab',false;
    'allSols',false;
    'msvOnly',true;
    'xplosRoots',false;
    'debug',false};

for ii=1:length(varargin)
    
    arg=varargin{ii};
    
    name=defaults{ii,1};
    
    if ~(islogical(arg) && isscalar(arg))
        
        error([name,' must be a logical and a scalar'])
    
    end
    
    defaults{ii,2}=arg;
    
end

slvOpts=cell2struct(defaults(:,2),defaults(:,1),1);

end
