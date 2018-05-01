function ab=optimize_logistic(x1,p1,x2,p2,op)
% OPTIMZE_LOGISTIC -- finds the optimal parameters of a logistic function
%
% ::
%
%
%   ab=optimize_logistic(x1,p1,x2,p2)
%
%   ab=optimize_logistic(x1,p1,x2,p2,op)
%
% Args:
%
%    - **x1,x2** [scalars] : points at which to evaluate the logistic
%
%    - **p1,p2** [scalars] : probability values for x1 and x2
%
%    - **op** ['+'|{'-'}] : decides whether b is to be multiplied by 1 or by
%      -1 in the evaluation the logistic function.
%
% Returns:
%    :
%
%    - **ab** [vector] : optimized parameters of the logistic function
%
% Note:
%
%      The function has the form f(x,a,b)=a/(a+exp(b*x)) so that f(x1,a,b)=p1
%      and f(x2,a,b)=p2
%
% Example:
%
%    See also:

if nargin<5
    
    op='-';
    
end
    
f=@(ab,x)ab(1)./(ab(1)+exp(ab(2)*x));

if strcmp(op,'-')
    
    op=@uminus;
    
elseif strcmp(op,'+')
    
    op=@(x)x;
    
else
    
    error('fifth argument should be + or -')
    
end

ab=fsolve(@objective,zeros(2,1),struct('Display','iter',...
    'TolFun',1e-12,'TolX',1e-07));

    function fval=objective(ab)
        
        ab(2)=op(ab(2));
        
        fval=[
            f(ab,x1)-p1
            f(ab,x2)-p2
            ];
        
    end

end