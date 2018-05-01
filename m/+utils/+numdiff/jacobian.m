function J=jacobian(func,x,varargin)
% jacobian - computes the jacobian of a function or a vector of functions
%
% ::
%
%
%   J=jacobian(func,x,varargin)
%
% Args:
%
%    - **func** [fhandle|cell array]: function or functions to be
%      differentiated
%
%    - **x** [n x 1 vector]: Vector of arguments for differentiation
%
%    - **varargin** : extra arguments of func beyond **x**
%
% Returns:
%    :
%
%    - **J** [matrix]: Numerical Jacobian of **func** at **x**
%
% Note:
%
%    The result of multiple functions is concatenated vertically.
%
% Example:
%
%    See also:

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


% try
%     
%     func{1}(xxp,varargin{:});
%     
%     col_by_col=false;
%     
% catch
%     
%     col_by_col=true;
%     
% end

J=cell(nrows,1);

for irow=1:nrows
    
%     if col_by_col
        
        for icol=1:n
            
            fi=func{irow}(xxp(:,icol),varargin{:})-...
                func{irow}(xxm(:,icol),varargin{:});
            
            if icol==1
                
                J{irow}=fi(:,ones(n,1));
                
            end
            
            J{irow}(:,icol)=fi;
            
        end
        
%     else
%         
%         J{irow}=func{irow}(xxp,varargin{:})-func{irow}(xxm,varargin{:});
%         
%     end
    
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

