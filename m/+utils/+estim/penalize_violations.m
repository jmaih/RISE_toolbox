function p=penalize_violations(viol,c,restrictions_same_weights)
% penalize_violations - penalizes the violation of general restrictions
%
% ::
%
%
%   p=penalize_violations(viol)
%
%   p=penalize_violations(viol,c)
%
%   p=penalize_violations(viol,c,restrictions_same_weights)
%
% Args:
%
%    - **viol** [vector] : value of constraints expressed in the form g(x)<=0
%
%    - **c** [scalar|{10}]: penalty factor
%
%    - **restrictions_same_weights** [true|{false}]: decides whether all
%    constraints are equivalent or not
%
% Returns:
%    :
%
%    - **p** [scalar] : penalty
%
% Note:
%
%    - restrictions are expected to be of the form g(x)<=0
%
% Example:
%
%    See also:

if nargin<3
    
    restrictions_same_weights=[];

if nargin<2
    
    c=[];
    
end

end

if isempty(restrictions_same_weights)
    
    restrictions_same_weights=false;
    
end

if isempty(c)
    
    c=10;
    
end

if c<=0
    
    error('penalty factor must be strictly positive')
    
end

p=0;

if ~isempty(viol)
    
    viol=max(0,viol(:));
    
    viol=viol(viol>0);
    
    if ~isempty(viol)
        
        if restrictions_same_weights
            
            p=numel(viol); % <-- sum(viol>0);
            
        else
        
        p=c*sum(viol.^2);% <--- p=c*sum((viol/min(viol)).^2);
        
        end
        
    end
    
end

end