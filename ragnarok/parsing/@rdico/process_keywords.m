%--- help for rdico.process_keywords ---
%
%  Without arguments, the function returns the list of model keywords
% 
%  Moving averages :: 
%  
%  	movave(X,n) = movavg(X,n)
%  	movavg(X,n) = (X{t} + X{t-1} + X{t-2} + ... + X{t-n+1})/abs(n)
%  	movave(X) = movavg(X) =	movave(X,-4)
%  
%  Moving products :: 
%  
%  	movprod(X,n) = prod(X{t},X{t-1},X{t-2},...X{t-n+1})
%  	movprod(X) = movprod(X,-4)
%  
%  Geometric averages :: 
%  
%  	movgeom(X,n) = prod(X{t},X{t-1},X{t-2},...X{t-n+1})^(1/abs(n))
%  	movgeom(X) = movgeom(X,-4)
%  
%  Moving sum:: 
%  
%  	movsum(X,n) = X{t} + X{t-1} + X{t-2} + ... + X{t-n+1}
%  	movsum(X) = movsum(X,-4)
%  
%  Differences:: 
%  
%  	diff(X,n) = X{t} - X{t+n}
%  	diff(X) = diff(X,-1)
%  
%  log differences ::
%  
%  	difflog(X,n) = log(X{t}) - log(X{t+n})
%  	difflog(X) = difflog(X-1)
%  
%  dot products ::
%  
%  	dot(X,n) = X{t}/X{t+n}
%  	dot(X) = dot(X,-1)
%