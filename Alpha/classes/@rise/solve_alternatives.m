function allobj=solve_alternatives(obj,varargin)
% syntax is allobj=solve_alternatives(obj,solver,file2save2)

if isempty(obj)
    if nargout>1
        error([mfilename,':: when the object is emtpy, nargout must be at most 1'])
    end
    allobj=struct('solve_alternatives_nsim',100,...
        'solve_alternatives_file2save2',[]);
    return
end

if ~isempty(varargin)
    obj=set_options(obj,varargin{:}); % <--- obj.options=mysetfield(obj.options,varargin{:});
end

solve_alternatives_nsim=obj.options.solve_alternatives_nsim;
file2save2=obj.options.solve_alternatives_file2save2;

% initialize the solution
allobj=rise.empty(0,0);
% load parameters from the mode, which can be empty of course
mode=obj.estimation.posterior_maximization.mode;
obj=assign_estimates(obj,mode);
[obj,retcode]=obj.solve();
sols=0;
indep_sols=0;
if ~retcode
    sols=sols+1;
    indep_sols=indep_sols+1;
    %     error([mfilename,':: model not solvable at the mode...'])
    allobj(1)=obj;
end
% for the rest, use random starts. In particular, we set the initial guess
% to true in order to start in different locations. This is not the ideal
% solution because I use rand, where all the random numbers are positive
% and in addition the magnitudes are not necessarily the ones in the "true"
% solution. Famer, Waggoner and Zha(2011) propose an algorithm for
% generating random starting values. I should try to adapt it...
obj.options.solve_initialization='random';
for sim=1:solve_alternatives_nsim-1
    obj=assign_estimates(obj,mode);
	[objtmp,retcode]=obj.solve();
    if ~retcode
        sols=sols+1;
        do_not_discard=true;
        for jj=1:indep_sols
            test=max(max(max(abs(objtmp.T-allobj(jj).T))));
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
obj.options.solve_initialization='default';
if ~isempty(file2save2)
    diary(file2save2)
end
for kk=1:indep_sols
    allobj(kk).print_solution
    if allobj(kk).is_stable_system
        disp('This system is stable')
    else
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