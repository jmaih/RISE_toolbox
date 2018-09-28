function [x,w] = quantization(n,xw_discrete)
% quantization -- discretizes a set of random variates
%
% ::
%
%   [x,w] = quantization(n,xw_discrete)
%
% Args:
%
%    n (vector): vector containing the number of discrete quantities in
%       each dimension of the gaussian shocks
%
%    xw_discrete (1 x 2 cell array): The first cell is the vector of all
%       possible values of the discrete distribution. The second cell is
%       the vector of weights of each of those possible values
%
% Returns:
%    :
%
%    - **x** [matrix]: grid combinations of the different variates, with
%       the number of shocks in rows and the number of combinations in
%       columns
%
%    - **w** [vector]: probability distribution over the different
%       combinations
%
% Note:
%
%    - Only works with gaussian shocks at the moment. 
%

if nargin<2
    
    xw_discrete=[];
    
end

is_discrete=~isempty(xw_discrete);

if is_discrete
    
    nd=numel(xw_discrete{1});
    
    n=[n(:).',nd];
    
end

d = length(n);

x = cell(1,d);

w0 = cell(1,d);

% the first elements are normal variates
for i=1:d-is_discrete
    
    [x{i},w0{i}] = quantize_normal(n(i));
    
end

% the last one is discrete
if is_discrete
    
    x{d}=xw_discrete{1}(:);
    
    w0{d}=xw_discrete{2}(:);
    
end

% invert in order to respect the kronecker structure of the grid below
w=utils.kronecker.kronall(w0{d:-1:1});

x = gridmake(x,n);

% transpose output
%-----------------
x=x.';

end

function [x,w] = quantize_normal(n)

% Based on an algorithm in W.H. Press, S.A. Teukolsky, W.T. Vetterling
% and B.P. Flannery, "Numerical Recipes in FORTRAN", 2nd ed.  Cambridge
% University Press, 1992.

maxit = 100;

pim4 = 1/pi.^0.25;

m = fix((n+1)./2);

x = zeros(n,1);

w = zeros(n,1);

for ii=1:m
    % starting values
    if ii==1
        
        z = sqrt(2*n+1)-1.85575*((2*n+1).^(-1/6));
        
    elseif ii==2
        
        z = z-1.14*(n.^0.426)./z;
        
    elseif ii==3
        
        z = 1.86*z+0.86*x(1);
        
    elseif ii==4
        
        z = 1.91*z+0.91*x(2);
        
    else
        
        z = 2*z+x(ii-2);
        
    end
    % root finding iterations
    iter=0;
    
    while iter<maxit
        
        iter = iter+1;
        
        p1 = pim4;
        
        p2 = 0;
        
        for j=1:n
            
            p3 = p2;
            
            p2 = p1;
            
            p1 = z.*sqrt(2/j).*p2-sqrt((j-1)/j).*p3;
            
        end
        
        pp = sqrt(2*n).*p2;
        
        z1 = z;
        
        z  = z1-p1./pp;
        
        if abs(z-z1)<1e-14
            
            break
            
        end
        
    end
    
    if iter>=maxit
        
        error('failure to converge in qnwnorm1')
        
    end
    
    x(n+1-ii) = z;
    
    x(ii) = -z;
    
    w(ii) = 2./(pp.*pp);
    
    w(n+1-ii) = w(ii);
    
end

w = w./sqrt(pi);

x = x*sqrt(2);

end

function x=gridmake(x0,n)

x=utils.gridfuncs.mygrid(n);

for ii=1:numel(n)
    
    x(:,ii)=x0{ii}(x(:,ii));
    
end

end
