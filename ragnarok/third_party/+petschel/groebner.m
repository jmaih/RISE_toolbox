function gbasis=groebner(polyset,ord,varnames,tol)
% GROEBNER - calculate the reduced Groebner basis of a set of polynomials
%  usage: gbasis = groebner(polyset,ord)
%         gbasis = groebner(polyset,ord,varnames)
%         gbasis = groebner(polyset,ord,varnames,tol)
%
% INPUTS: polyset, ord, varnames, tol
%  polyset is a cell array of multivariate polynomial coefficients
%    the cell elements can be all strings or all arrays.
%    If the polynomials are specified with arrays, each row of the array
%    represents a single term.  The first element is the coefficient and
%    the others are the degrees of each unknown in the term.
%    e.g. [-1,2,1;1,0,1] represents -x1^2*x2+x2
%    If the polynomials are specified as strings, they must use variable
%    names 'x1','x2',... (unless varnames is provided) and be valid inputs
%    for petschel.str2poly.m
%  ord is a string specifying the monomial ordering:
%    'lex': lexicographical, order by highest power of most significant
%      indeterminate, e.g. x1>x2^2, x1*x2>x1.
%    'grlex': graded lex, order by total degree then lex,
%      e.g. x1^2*x2 > x1*x2^2 > x1^2
%    'grevlex': graded reverse lex, order by total degree then by lowest
%      power of least significant indeterminate, e.g. x1*x2^2 > x1*x2*x3
%    (for all cases have x1>x2>...>xn)
%    note: revlex is not a well-ordering because 1>x1>x1^2>...
%  varnames (optional) is a cell array of variable names if polyset is a
%    cell array of strings and the variable names are not 'x1', 'x2', etc.
%    If specified, there must always be at least 2 variables.
%  tol (optional) is the zero tolerance for coefficients, otherwise
%    catastrophic cancellation may occur.  Default value is 0
%
% OUTPUTS: gbasis
%  gbasis is the reduced Groebner basis of the set of polynomials, in the
%    same format (strings or matrices) as polyset{1}.
%
% ALGORITHM:
%  A modified Buchberger's algorithm is used.  In Buchberger's algorithm,
%  at each step calculate Sij = (Lij/ai)*Pi - (Lij/aj)*Pj, for all pairs of
%  polynomials Pi, Pj, where ai and aj are the leading terms of Pi and Pj,
%  and Lij is the least common multiple of ai and aj; then reduce Sij with
%  respect to the polynomials and add it to the set (note that Sij always
%  reduces to zero if ai and aj have no variables in common).  The
%  algorithm terminates when no new polynomials can be added to the set.
%  This function reduces the polynomials at each step in order to find the
%  reduced Groebner basis.
%
% EXAMPLE: simplify equations x^2+2xy^2=0, xy+2y^3=1
%  A1=[1,2,0;2,1,2]; % 1*x^2*y^0 + 2*x^1*y^2
%  A2=[-1,0,0;1,1,1;2,0,3]; % -1*x^0*y^0 + 1*x^1*y^1 + 2*x^0*y^3
%  y=groebner({A1,A2},'lex');
%  petschel.poly2str(y) % returns {'x2^3-0.5', 'x1'}
%
% EXAMPLE: same, but with equations specified as strings:
%  groebner({'x^2+2*x*y^2','x*y+2*y^3-1'},'lex',{'x','y'})
%  % returns {'y^3-0.5','x'}
%
% EXAMPLE: show that a set of equations is inconsistent
%  groebner({'x + y', 'x^2 - 1', 'y^2 - 2*x'}, 'lex', {'x', 'y'})
%  % returns {'1'}
%
%
% KNOWN ISSUES
%
%  0. GENERAL: while a significant effort has been made to ensure that the
%    program works as described, unidentified bugs may remain and the results
%    should be regarded with caution.  In particular, the suite of tests
%    performed was by no means comprehensive.
%
%  1. ACCURACY: the algorithm by default works in floating-point which
%    ultimately cannot perform exact arithmetic, so unexpected silent errors may
%    occur even when the coefficients are specified as integers.  A workaround
%    is to use other data types but this has its own issues (see below).
%
%  1a. EXAMPLE (roundoff - linearly independent polynomials)
%    eps % returns 2.2204e-16 (32-bit PC windows machine epsilon)
%    groebner({'x1+x2','x1+1.000000000000001*x2'},'lex') % ok, returns {'x2','x1'}
%    groebner({'x1+x2','x1+1.0000000000000001*x2'},'lex') % not ok, returns {'x2+x1'}
%
%  1b. EXAMPLE (roundoff - linearly dependent polynomials)
%    2/sqrt(2)-sqrt(2) % expect 0, returns -2.2204e-016
%    petschel.poly2str(groebner({[1,1,0;sqrt(2),0,1],[1/sqrt(2),1,0;1,0,1]},'lex')) % ok, returns {'1.41421*x2+x1'}
%    petschel.poly2str(groebner({[1,1,0;sqrt(2),0,1],[sqrt(2),1,0;2,0,1]},'lex')) % not ok, returns {'x2','x1'}
%    petschel.poly2str(groebner({[1,1,0;sqrt(2),0,1],[sqrt(2),1,0;2,0,1]},'lex',{},1e-14)) % ok, returns {'1.41421*x2+x1'}
%
%  1c. EXAMPLE (courtesy of Christophe Lauwerys - catastrophic cancellation & tolerance bug)
%    groebner({'t^3+x+y','t^2+0.5*x^2-x-z^2','t^2+y-z^2'},'lex',{'t','x','y','z'}) % ans{1} contains Inf
%
%  1d. EXAMPLE (from Mathematica documentation - catastrophic cancellation because of large intermediate coefficients)
%    % expect 3 polynomials returns, one of degree 8 in z; instead can get '1'
%    groebner({'x^2 + y^2 + z^2 - 1', 'x*y - z + 2', 'z^2 - 2*x + 3*y'}, 'lex', {'x', 'y','z'})
%
%  2. NON-DOUBLE INPUTS: use of non-double coefficients is not supported.  While
%    it may be tempting to introduce symbolic coefficients, it is not clear how
%    they should be consistently treated - for example 'A*x' could reduce to '0'
%    or 'x' or both, depending upon the assumptions on the coefficient A.
%
%  3. EFFICIENCY: The reduction algorithms aren't terribly efficient.  One of
%    the major bottlenecks at present is the calculation of the leading term,
%    which currently involves a call to sortrows each time.  It could be more
%    efficient to maintain an ordered list or even a priority queue (e.g.
%    implemented as a heap).  Also lattice reduction algorithms such as LLL
%    might speed things up.
%
%
% SEE ALSO:
%  petschel.poly2str, petschel.str2poly, petschel.polynsolve

% Author: Ben Petschel 19/6/2009
%
% Change history:
%  19/6/2009 - first release (ord='grlex' and ord='grevlex' not
%    implemented yet)
%  22/6/2009 - implemented ord='grlex' and ord='grevlex'
%    reduction algorithm now reduces lower-order terms as well
%  23/6/2009 - removed ord='revlex' because it is not a well-ordering
%  17/7/2009 - fixed bug in handling polynomials with >=3 variables
%  20/3/2010 - changed polynomial representation from multidimensional
%    arrays to 2d arrays (more memory efficient with large number of vars)
%  11/10/2010 - bugfix for handling error tolerance

if (numel(polyset)>0) && ischar(polyset{1})
    
    charout = true;
    
    if nargin<3
        
        polyset = petschel.str2poly(polyset);
        
    else
        
        polyset = petschel.str2poly(polyset,varnames);
        
    end
    
else
    
    charout = false;
    
end

if nargin<4
    
    tol = 0;
    
end

% determine the number of degrees
nd=max(cellfun(@(x)size(x,2)-1,polyset)); % error if = 0

if nd == 0
    
    error('number of variables must be >= 1');
    
end

gbasis = fullreduce(polyset,ord,tol,nd);

oldgbasis = {};

while ~isequal(oldgbasis,gbasis)
    
    oldgbasis = gbasis;
    
    Qset = SPoly(gbasis,ord,tol,nd);
    
    gbasis = fullreduce([gbasis,Qset],ord,tol,nd);
    
end

if charout
    
    if nargin<3
        
        gbasis = petschel.poly2str(gbasis);
        
    else
        
        gbasis = petschel.poly2str(gbasis,varnames);
        
    end
    
end

end % main function groebner(...)


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function Qset = SPoly(Pset,ord,tol,nd)
% calculates all reduced S-polynomials of Pset (Pset must be reduced)
%   Sij = (Lij/ai)*Pi - (Lij/aj)*Pj,
% for all pairs of polynomials Pi, Pj, where ai and aj are the leading
% terms of Pi and Pj, and Lij is the least common multiple of ai and aj
% then reduces Sij wrt Pset

Qset = {};

for i=1:numel(Pset)-1
    
    for j=i+1:numel(Pset)
        
        [c1,d1] = leadterm(Pset{i},ord,tol,nd); %#ok<*ASGLU>
        
        [c2,d2] = leadterm(Pset{j},ord,tol,nd);
        
        if any(and(d1>0,d2>0))
            % leading terms have a variable in common, so S-poly is nontrivial
            L = max(d1,d2); % LCM of leading terms
            
            S = addpolys(multiplyterm(Pset{i},L-d1),multconst(multiplyterm(Pset{j},L-d2),-1));
            
            S = reduceset(S,Pset,ord,tol,nd);
            
            %if any(S(:)~=0),
            if abs(leadterm(S,ord,tol,nd))>tol
                
                Qset{end+1}=S; % add S to the list
                
            end
            
        end
        
    end
    
end

end % helper function SPoly(...)


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function Qset = fullreduce(Pset,ord,tol,nd)
% reduces all polynomials wrt each other
Qset = Pset;

oldQset = {};

while ~isequal(oldQset,Qset)
    
    oldQset = Qset;
    
    for i=1:numel(Qset)-1
        
        for j=i+1:numel(Qset)
            
            Qset{i}=reduce(Qset{i},Qset{j},ord,tol,nd);
            
            Qset{j}=reduce(Qset{j},Qset{i},ord,tol,nd);
            
        end
        
    end
    
end

% keep only the nonzero results (to within zero tolerance tol)
Qset = Qset(cellfun(@(x)abs(leadterm(x,ord,tol,nd))>tol,Qset));

if numel(Qset)==1
    
    % special case, a single equation can slip through without being reduced
    c = leadterm(Qset{1},ord,tol,nd);
    
    % coefficients occupy first column of the array
    Qset{1} = multconst(Qset{1},1/c); % know c~=0 because have already removed zero eqns
    
end

end % helper function fullreduce(...)


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function Q = reduceset(P,Pset,ord,tol,nd)
% reduces P wrt all polynomials in Pset
Q = P;

oldQ = [];

while ~isequal(oldQ,Q)
    
    oldQ = Q;
    
    for i=1:numel(Pset)
        
        Q=reduce(Q,Pset{i},ord,tol,nd);
        
    end
    
end

end % helper function reduceset(...)


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function Q1=reduce(P1,P2,ord,tol,nd)
% reduces a polynomial P1 by subtracting multiples of P2 so that leading
% term coefficient is 1 and leading term of P2 does not divide that of P1

Q2=P2;

[c2,d2] = leadterm(Q2,ord,tol,nd);

if abs(c2)<=tol
    
    % ignore terms less than tol
    keepgoing = false;
    
    Q1 = P1;
    
else
    
    keepgoing = true;
    
    Q2 = multconst(Q2,1/c2);
    
    Q1 = 0; % will successively add terms to Q
    % remove the leading terms, in order to reduce wrt lower terms
    remain = P1;
    
end


while keepgoing
    
    [c1,d1] = leadterm(remain,ord,tol,nd);
    
    if abs(c1)<=tol
        
        % ignore terms less than tol
        keepgoing = false;
        
    else
        %reduce remain by Q2, if possible, or remove leading term
        if d1>=d2
            % can reduce Q1 by Q2, so subtract a multiple of Q2 from Q1
            remain = addpolys(remain,multconst(multiplyterm(Q2,d1-d2),-c1));
            
        else
            % cannot reduce any further wrt leading term, but try lower terms
            lead = multiplyterm(c1,d1);
            
            remain = addpolys(remain,multconst(lead,-1));
            
            Q1 = addpolys(Q1,lead); % add irreducible terms to Q1
            
        end
        
    end
    
end

% reduce leading coefficient, if possible
c1 = leadterm(Q1,ord,tol,nd);

if abs(c1)>tol
    
    Q1 = multconst(Q1,1/c1);
    
end

end % helper function reduce(...)


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function Q=multconst(P,c)
% returns c*P where c is a constant

Q = P;

Q(:,1)=Q(:,1)*c;

end % helper function multconst(...)


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function Q=addpolys(P1,P2)
% adds two polynomials

s1 = size(P1,2);

s2 = size(P2,2);

if s1>s2
    
    Q = [P1;P2,zeros(size(P2,1),s1-s2)];
    
else
    
    Q = [P1,zeros(size(P1,1),s2-s1);P2];
    
end

Q = sortrows(Q,2:size(Q,2));

i=1;

while i<size(Q,1)
    % merge rows representing terms with the same exponents
    if Q(i,1)==0
        
        Q(i,:)=[];
        
    elseif all(Q(i,2:end)==Q(i+1,2:end))
        
        Q(i,1)=Q(i,1)+Q(i+1,1);
        
        Q(i+1,:)=[];
        
        if Q(i,1)==0
            
            % remove terms that add to zero
            Q(i,:)=[];
            
        end
        
    else
        
        i=i+1;
        
    end
    
end


end % helper function addpolys


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function Q=multiplyterm(P,d)
% multiplies a polynomial by a term (no coefficient)

% padding a matrix at the end is same as multiplying by x
d=d(:).';

sP = size(P,2);

sd = size(d,2)+1; % d doesn't include coefficient

if sP<sd
    
    Q=[P,zeros(size(P,1),sd-sP)];
    
else
    
    Q=P;
    
    d=[(d(:)).',zeros(1,sP-sd)];
    
end

Q=Q+repmat([0,d],size(Q,1),1);

end % helper function multiplyterm


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [coeff,deg]=leadterm(P,ord,tol,nd)
% returns the leading term of P
% deg is a row vector of length size(P,2)-1
% coeff is the leading coefficient, deg(i) is the degree of xi
% ignores terms with coefficients <= tol
% nd is the minimum number of degrees (in order to handle cases where deg=0)

d=size(P,2);

if all(abs(P(:,1))<=tol)
    % P=0
    coeff=0;
    deg=zeros(1,d-1);
else
    switch ord
        case 'lex'
            P=sortrows(P,2:d);
            
        case 'revlex'
            error('revlex is not a well-ordering');
            
        case 'grlex'
            P = [P,sum(P(:,2:d),2)];
            P = sortrows(P,[d+1,2:d]);
            P = P(:,1:d);
            
        case 'grevlex'
            P = [P,sum(P(:,2:d),2)];
            P = sortrows(P,[d+1,-(d:-1:2)]);
            P = P(:,1:d);
            
        otherwise
            error('ordering %s not recognized',ord);
            
    end % switch ord
    
    ind = find(abs(P(:,1))>tol,1,'last');
    
    if isempty(ind)
        
        coeff = 0;
        
        deg = zeros(1,d-1);
        
    else
        
        coeff = P(ind,1);
        
        deg = P(ind,2:end);
        
    end
    
end % if all(P(:)==0)

if nd>d-1
    
    deg = [deg,zeros(1,nd-d+1)];
    
end

end % helper function leadterm(...)


