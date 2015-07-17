function J=jacobian(func,x,varargin)
% jacobian - computes the jacobian of a function or a vector of functions
%
% Syntax
% -------
% ::
%
%   J=jacobian(func,x,varargin)
%
% Inputs
% -------
%
% - **func** [fhandle|cell array]: function or functions to be
%   differentiated
%
% - **x** [n x 1 vector]: Vector of arguments for differentiation
%
% - **varargin** : extra arguments of func beyond **x**
%
% Outputs
% --------
%
% - **J** [matrix]: Numerical Jacobian of **func** at **x**
%
% More About
% ------------
%
% The result of multiple functions is concatenated vertically.
%
% Examples
% ---------
%
% See also: 

tol=eps^(1/3);

x=x(:);
h=tol*max(1,x);
xp=x+h;
xm=x-h;
h=xp-xm;
n=numel(x);

x=x(:,ones(n,1));

diag_terms=(0:n-1)*n+(1:n);
xxp=x;
xxp(diag_terms)=xp;
xxm=x;
xxm(diag_terms)=xm;

if ~iscell(func)
    func={func};
end
[nrows,ncols]=size(func);
if ncols>1
    error('the objective cannot have many columns since outputs are to be concatenated vertically')
end

J=cell(nrows,1);
for irow=1:nrows
    J{irow}=func{irow}(xxp,varargin{:})-func{irow}(xxm,varargin{:});
    nrows_i=size(J{irow},1);
    hp=h';
    J{irow}=J{irow}./hp(ones(nrows_i,1),:);
end

if nrows==1
    J=J{1};
else
    J=cell2mat(J);
end

end

