function block=cell2block(tmp,pname,chain,state,name_file_line)

n_entries=numel(tmp);

if ~iscell(tmp)|| n_entries < 3
    
    error('all fields of the prior structure should be cell arrays with at least 3 elements')
    
end

start=parser.push_if_validated(tmp{1},@(x)isfinite(x),'start value',name_file_line);

lq=nan;    uq=nan;     pmean=nan;    pstdev=nan;    prior_prob=1;

distrib='uniform';    lb=nan;    ub=nan;

if n_entries==3
    
    lb=parser.push_if_validated(tmp{2},@(x)isfinite(x),'lower bound',name_file_line);
    
    ub=parser.push_if_validated(tmp{3},@(x)isfinite(x)&& x>lb,'upper bound value',name_file_line);
    
    lq=lb;
    
    uq=ub;
    
else
    
    distrib=parser.push_if_validated(tmp{4},@(x)ischar(x),'distribution(prob)',name_file_line);
    
    left_par=strfind(distrib,'(');
    
    if isempty(left_par)
        
        left_par=length(distrib)+1;
        pmean=parser.push_if_validated(tmp{2},@(x)isfinite(x),'prior mean',name_file_line);
        
        pstdev=parser.push_if_validated(tmp{3},@(x)isfinite(x)&& x>0,'prior stdev',name_file_line);
        
        % we have to change default for the probability
        prior_prob=nan;
        
    else
        right_par=strfind(distrib,')');
        
        prior_prob=eval(distrib(left_par+1:right_par-1));
        
        prior_prob=parser.push_if_validated(prior_prob,@(x)isfinite(x) && x>=0 && x<=1,'prior probability',name_file_line);
        
        lq=parser.push_if_validated(tmp{2},@(x)isfinite(x),'lower quantile',name_file_line);
        
        uq=parser.push_if_validated(tmp{3},@(x)isfinite(x)&& x>lq,'upper quantile',name_file_line);
        
    end
    
    distrib=distrib(1:left_par-1);
    
    if n_entries>4
        
        lb=parser.push_if_validated(tmp{5},@(x)isfinite(x),'lower bound',name_file_line);
        
        if n_entries>5
            
            ub=parser.push_if_validated(tmp{6},@(x)isfinite(x),'upper bound',name_file_line);
            
            if n_entries>6
                
                error('number of entries in setting up the prior cannot exceed 6')
                
            end
            
        end
        
    end
    
end

block=struct('name',pname,'chain',chain,'state',state,'start',start,...
    'lower_quantile',lq,'upper_quantile',uq,'prior_mean',pmean,...
    'prior_stdev',pstdev,'prior_distrib',distrib,'prior_prob',prior_prob,...
    'lower_bound',lb,'upper_bound',ub);

end