function [x,f,retcode]=generate_starting_point(objective,max_trials)
% H1 line
%
% ::
%
%
% Args:
%
% Returns:
%    :
%
% Note:
%
% Example:
%
%    See also:

% this function attempts to reduce the number of rejection of randomly
% drawn vectors in the presence of linear or nonlinear restrictions. When
% there are many restrictions, randomly drawing a vector of parameters that
% satisfies them all can be challenging. The strategy is to select the
% vector with the least possible constraints.
if nargin<2
    
    max_trials=[];
    
end

if isempty(max_trials)
    
    max_trials=50;
    
end

[x,f,~,retcode,viol]=objective();

iter=0;

while iter<max_trials && sum(viol)>0
    
    iter=iter+1;
    
    [newc,fc,~,retcodec,violc]=objective();
    
    if sum(violc)<sum(viol)
        
        x=newc;
        
        viol=violc;
        
        retcode=retcodec;
        
        f=fc;
        
    end
    
end

end