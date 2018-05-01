function [allobj,mean_square_stable]=solve_alternatives(obj,varargin)
% solve_alternatives - attempts to find all solutions of a MSDSGE model
%
% ::
%
%   [allobj,mean_square_stable]=solve_alternatives(obj)
%   [allobj,mean_square_stable]=solve_alternatives(obj,varargin)
%
% Args:
%
%    - **obj** [rise|dsge]: RISE model object
%
%    - **solve_alternatives_nsim** [numeric|{100}]: Number of random starting
%    values generated.
%
%    - **solve_alternatives_file2save2** [char|{''}]: name of the diary file
%    where the results are saved.
%
% Returns:
%    :
%
%    - **allobj** [rise|dsge]: vector of RISE objects with the solutions found
%
%    - **mean_square_stable** [true|false]: vector of logicals with the same
%    number of elements as allobj. An entry equal to true means the particular
%    model is mean square stable.
%
% Note:
%
% Example:
%
%    See also:

% syntax is allobj=solve_alternatives(obj,solver,file2save2)

if isempty(obj)
    
    mydefaults=the_defaults();
    
    if nargout
        
        allobj=mydefaults;
        
    else
        
        disp_defaults(mydefaults);
        
    end
    
    return
    
end

if ~isempty(varargin)
    
    obj=set(obj,varargin{:}); 
    
end

solve_alternatives_nsim=obj.options.solve_alternatives_nsim;

file2save2=obj.options.solve_alternatives_file2save2;

% initialize the solution
%-------------------------
classobj=class(obj);

allobj=eval([classobj,'.empty(0)']);

[obj,retcode]=obj.solve();

sols=0;

indep_sols=0;

if ~retcode
    
    sols=sols+1;
    
    indep_sols=indep_sols+1;
    
    allobj(1)=obj;
    
end

iter=sols;

obj=set(obj,'solve_initialization','random');

for sim=iter+1:solve_alternatives_nsim
    
	[objtmp,retcode]=obj.solve();
    
    if ~retcode
        
        sols=sols+1;
        
        do_not_discard=true;
        
		newsol=cell2mat(objtmp.solution.Tz);
        
        for jj=1:indep_sols
            
            test=max(max(abs(newsol-cell2mat(allobj(jj).solution.Tz))));
            
            if test<1e-5
                
                do_not_discard=false;
                
                break
                
            end
            
        end
        
        if do_not_discard
            
            indep_sols=indep_sols+1;
            
            allobj(indep_sols,1)=objtmp;
            
        end
        
    end
    
end

if ~isempty(file2save2)
    
    diary(file2save2)
    
end

mean_square_stable=true(1,indep_sols);

for kk=1:indep_sols
    
    allobj(kk).print_solution
    
    if allobj(kk).is_stable_system
        
        disp('This system is stable')
        
    else
        
        mean_square_stable(kk)=false;
        
        disp('This system is not stable')
        
    end
    
end

disp([mfilename,':: for ',int2str(solve_alternatives_nsim),' starting points, ',...
    int2str(sols),' solutions were found, but only ',...
    int2str(indep_sols),' were independent'])

if ~isempty(file2save2)
    
    diary off
    
end

end

function d=the_defaults()

num_fin=@(x)isnumeric(x) && isscalar(x) && isfinite(x);

num_fin_int=@(x)num_fin(x) && floor(x)==ceil(x) && x>=0;

d={
    'solve_alternatives_nsim',100,@(x)num_fin_int(x),...
    'solve_alternatives_nsim must be a finite and positive integer'
    
    'solve_alternatives_file2save2',[],@(x)ischar(x),...
    'solve_alternatives_file2save2 must be a char'
    };
    
end