function V = multivariate_taylor_approximation(ghandle,values,order,varargin)
%MVT main program: returns a cell array of all multivariate Taylor Coefficients
%   of g given by @g_func_name, about values (1xn), through degree order.
%   by Richard Neidinger, 6/4/04.  This is Algorithm 1 of his paper
%  "Directions for Computing Truncated Multivariate Taylor Series"

n = length(values);  % number of variables
[W,Binom,B]=makeGlobals(n,order);

% initialize x to be a vector of series objects for each variable.
x = dyn_series.empty(0,0);
for i = 1:n
    x(1,i)=dyn_series(values(i),1,order);
end

V = cell(order+1,1);  % V{k} will be all (k-1)th order series terms from all needed directions, then,
% divDiff and collect_terms will interpolate to (k-1)th order multivariate coefficients.
for k = 1:order+1
    alpha = [k-1 zeros(1,n-2)];
    while alpha(1)<k
        x = assigncoef1(x,[1 W(alpha+1)]);  % changes direction in variables.
        uvec = double(ghandle(x,varargin{:}));  % automatic differentiation using @series
        for j = k:order+1
            V{j} = [V{j},uvec(:,j)];
        end
        alpha = increment(alpha);
    end
    V{k} = divDiff(V{k},n-1,k-1,W,Binom);
    V{k} = collect_terms(V{k},n-1,k-1,Binom,B);
end

function [W,Binom,B]=makeGlobals(nmax,dmax)
% MAKEGLOBALS  creates global arrays W, Binom, and B
%   W is used by for use by mvt.m, divDiff.m, collect_terms.m;
%      to access w_i in the paper use W(i+1) since starting index is 1, not 0.
%   Binom is used by index.m;  Binom(i,j) = Pascal[j,i-2] from paper.
%   B is used by collect_terms.m;  B(i,j) = b_(i,j) in paper.


% alternating dyadic rationals for W node values
W = [0 1 -1];
i = 1;
while length(W)<dmax+1
    p = 2^i;
    odds = 1:2:(p-1);
    W = [W,(1/p)*reshape([odds;-odds],1,p)];
    i = i+1;
end
%  W = 0:(dmax+1);   % alternative of nonnegative integers for W node values

%  Binom(i,j) = ((i+j-2) choose j) = Pascal[j,i-2] from paper
Binom = (0:dmax)';
for i = 1:nmax-1
    Binom = [Binom, cumsum(Binom(:,end))];
end

%  B array from paper for i,j > 0
%  Assumes that w_0 = W(1) is zero!
B = ones(dmax,dmax);
B(1,(2:dmax)) = cumprod(-W(2:dmax));
for j = 3:dmax
    B(2:(j-1),j) = B(1:(j-2),j-1) - W(j)*B(2:(j-1),j-1);
end

function V = divDiff(V,n,k,W,Binom)
%divDiff takes a vector of wedge-node func-values for n variables thru degree k.
%   Returns vector of coefficients of interpolating newton polynomial.
%   Uses globals W, whose index starts with one, not zero as
%   by Richard Neidinger, 6/7/04.  This is Algorithm 3 of his paper
%  "Directions for Computing Truncated Multivariate Taylor Series"

% W  node w0,w1,w2,...,wd sequence defined with global W before running this.

endi = size(V,2);
if endi ~= index([zeros(1,n-1), k],Binom)
    disp('Error in input dimensions for divDiff.');
    return;
end
for jj=1:k % pass j through the tableau
    beta = [zeros(1,n-1), k];
    i = endi;
    jstarti = index([jj, zeros(1,n-1)],Binom);
    while i >= jstarti
        m = 1;
        psum = beta(1);
        while psum < jj
            m = m+1;
            psum = psum + beta(m);
        end  % psum
        b = beta(m);
        a = psum - jj;
        bdown = beta;
        bdown(m) = b-1;
        V(:,i) = (V(:,i) - V(:,index(bdown,Binom))) / (W(b+1) - W(a+1));
        beta = decrement(beta);
        i = i-1;
    end
end  % j pass

function isum = index(alpha,Binom)
%INDEX takes a multi-index vector and returns the scalar index for it
%   within a list, starting with one, of all multi-indices in order.
%   by Richard Neidinger, 6/8/04.  This is (7.4) of his paper
%  "Directions for Computing Truncated Multivariate Taylor Series"

n = length(alpha);
p=1;
isum = 1;
for j = 1:n
    p = p + alpha(n+1-j);
    isum = isum + Binom(p,j);
end

function V = collect_terms(V,n,k,Binom,B)
%collect_terms takes a vector of coefficients of interpolating newton polynomial
%   for n variables thru degree k.  (Paper uses n-1 in place of n here.)
%   Returns vector of coefficients of standard expanded polynomial form.
%   by Richard Neidinger, 6/9/04.  This is Algorithm 5 of his paper
%  "Directions for Computing Truncated Multivariate Taylor Series"

endi = size(V,2);
if endi ~= index([zeros(1,n-1), k],Binom)
    disp('Error in input dimensions for collect_terms.');
    return;
end
for order = 1:k-1
    alpha = [order, zeros(1,n-1)];
    i = index(alpha,Binom);
    while alpha(1) <= order
        nz = find(sparse(alpha));  % vector of locations of nonzero elements
        m = length(nz);
        mu = [1, zeros(1,m-1)];
        while mu(1) <= k - order
            beta = alpha;
            beta(nz) = beta(nz) + mu;
            bprod = 1;
            for c = nz
                bprod = bprod * B(alpha(c),beta(c));
            end
            V(:,i) = V(:,i) + bprod * V(:,index(beta,Binom));
            mu = increment(mu);
        end  % mu
        alpha = increment(alpha);
        i = i+1;
    end  % alpha
end  % order

function alpha = decrement(alpha)
%DECREMENT takes a multi-index vector and returns the previous multi-index.
%   input should be a vector of nonnegative integers.
%   by Richard Neidinger, 6/8/04.
m = length(alpha);
i = m;
while alpha(i)==0
    i = i-1;
end
alpha(m) = alpha(i) - 1;
if i < m
    alpha(i) = 0;
end
if i>1
    alpha(i-1) = alpha(i-1) + 1;
end

function alpha = increment(alpha)
%INCREMENT takes a multi-index vector and returns the next multi-index.
%   input should be a vector of nonnegative integers.
%   by Richard Neidinger, 6/3/04.  This is Algorithm 2 of his paper
%  "Directions for Computing Truncated Multivariate Taylor Series"
m = length(alpha);
i = m-1;
while i>0 && alpha(i)==0
    i = i-1;
end
if i>0
    alpha(i) = alpha(i)-1;
end
alpha(i+1) = alpha(m) + 1;
for j = (i+2):m
    alpha(j) = 0;
end