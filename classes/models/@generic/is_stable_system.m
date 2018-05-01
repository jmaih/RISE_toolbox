function flag=is_stable_system(obj,varargin)
% IS_STABLE_SYSTEM -- checks the stability of a linear markov switching
% system.
%
% ::
%
%
%   flag=IS_STABLE_SYSTEM(obj)
%
%   flag=IS_STABLE_SYSTEM(obj,varargin)
%
% Args:
%
%    - **obj** [dsge|rise|svar|rfvar]: model object
%
%    - **varargin** [name,value]: pairwise valid options for RISE. The most
%    relevant in this case are
%      - **stability_criterion** [numeric|{1.000001}]: stability criterion.
%      All eigenvalues must be smaller than this criterion for the system to
%      be MSS
%      - **stability_algorithm** ['cfm'|{'hmg'}]: CFM stands for
%      Costa-Fragoso-Marques while HMG stands for Hassibi-Murray-Gupta.
%
% Returns:
%    :
%
%    - **flag** [false|true]: result of the investigation on whether the
%    system is stable or not.
%
% Note:
%
% Example:
%
%    See also:

% this function checks that the solved system is
% stable in the sense that its covariance matrix is
% bounded.
if isempty(obj)
    
    mydefaults=the_defaults();
    
    if nargout
        
        flag=mydefaults;
        
    else
        
        disp_defaults(mydefaults);
        
    end
    
    return
    
end

if ~isempty(varargin)
    
    obj=set(obj,varargin{:});

end

nobj=numel(obj);

if nobj>1
    
    flag=nan(1,nobj);
    
    for iobj=1:nobj
        
        flag(iobj)=is_stable_system(obj(iobj));
    
    end
    
    return

end
% this is hard-coded for the moment. We don't want to check stability under
% estimation if the problem is too big...
ms_stability_check_threshold=5000;

% trim down in case some matrices are the same and keep only the states
[T,Q,n,h]=problem_reduction(obj);

flag=true;

if ~isempty(T{1})
    % update n right here right now
    n2=n^2;
    
    if ~(obj.estimation_under_way && h*n2>ms_stability_check_threshold)
        
        crit=obj.options.stability_criterion;
        % do not check stability under estimation if matrix is too big. I
        % should probably decrease the threshold coz what is expensive is
        % the calculation of the kronecker products. Moreover, I could put
        % the threshold as an option...
        switch lower(obj.options.stability_algorithm)
            
            case 'gmh'
                
                flag=utils.mss.gupta_murray_hassibi(T,Q,crit);
            
            case 'cfm'
                
                flag=utils.mss.costa_fragoso_marques(T,Q,crit);
            
            otherwise
                
                error('valid stability algorithms are gmh and cfm')
        
        end
        
    end
    
end

end

function d=the_defaults()

num_fin=@(x)isnumeric(x) && isscalar(x) && isfinite(x);

d={
    'stability_criterion',1.000001,@(x)num_fin(x) && x>=1,...
    'stability_criterion must be >=1'
    
    'stability_algorithm','gmh',@(x)ismember(x,{'gmh','cfm'}),...
    'stability_algorithm must be ''gmh'' or ''cfm'''
    };

end