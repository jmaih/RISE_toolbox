function allobj=solve_alternatives(obj,varargin)
% Attempts to find all solutions of a MSDSGE model
%
% ::
%
%   [allobj,mean_square_stable] = solve_alternatives(obj)
%   [allobj,mean_square_stable] = solve_alternatives(obj,varargin)
%
% Args:
%
%    obj (rise | dsge): RISE model object
%    solve_alternatives_nsim (numeric | {100}): Number of random starting
%      values generated.
%
% Returns:
%    :
%
%    - **allobj** [rise|dsge]: vector of RISE objects with the solutions found

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

% initialize the solution
%-------------------------
obj.solution=[];

obj=set(obj,'solve_check_stability',false,...
    'solve_initialization','random');

sols=0;

indep_sols=0;

tol=1e-5;

for sim=1:solve_alternatives_nsim
    
	[objtmp,retcode]=solve(obj);
    
    if ~retcode
        
        if sols==0
            
            solution_map=objtmp.solution;
            
        end
        
        sols=sols+1;
        
        do_not_discard=true;
        
        newmap=objtmp.solution;
        
		newsol=cell2mat(newmap.Tz);
        
        for jj=1:indep_sols
            
            test=max(max(abs(newsol-cell2mat(solution_map(jj).Tz))));
            
            if test<tol
                
                do_not_discard=false;
                
                break
                
            end
            
        end
        
        if do_not_discard
            
            indep_sols=indep_sols+1;
            
            solution_map(indep_sols)=newmap;
            
        end
        
    end
    
end

allobj=obj;

if indep_sols==0
    
    return
    
end
    
finalSol=solution_map(1);

zzz=repmat('z',1,obj.options.solve_order);

for kk=2:indep_sols
    
    for io=1:obj.options.solve_order
        
        tzz=['T',zzz(1:io)];
        
        finalSol.(tzz)=cat(3,finalSol.(tzz),solution_map(kk).(tzz));
        
    end
    
end

allobj.solution=finalSol;

disp([mfilename,':: for ',int2str(solve_alternatives_nsim),' starting points, ',...
    int2str(sols),' solutions were found, but only ',...
    int2str(indep_sols),' were independent'])

end

function d=the_defaults()

num_fin=@(x)isnumeric(x) && isscalar(x) && isfinite(x);

num_fin_int=@(x)num_fin(x) && floor(x)==ceil(x) && x>=0;

d={
    'solve_alternatives_nsim',100,@(x)num_fin_int(x),...
    'solve_alternatives_nsim must be a finite and positive integer'
    };
    
end