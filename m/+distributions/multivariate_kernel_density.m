function [F,X]=multivariate_kernel_density(X0,N,kernel,id,findex,xindex,F,X,Hi_dHi)
% multivariate_kernel_density -- multivariate kernel density estimation
%
% Syntax
% -------
% ::
%
%   [F,X]=multivariate_kernel_density(X0)
%
%   [F,X]=multivariate_kernel_density(X0,N)
%
%   [F,X]=multivariate_kernel_density(X0,N,kernel)
%
%   [F,X]=multivariate_kernel_density(X0,N,kernel,id)
%
%   [F,X]=multivariate_kernel_density(X0,N,kernel,id,findex)
%
%   [F,X]=multivariate_kernel_density(X0,N,kernel,id,findex,xindex)
%
%   [F,X]=multivariate_kernel_density(X0,N,kernel,id,findex,xindex,F)
%
%   [F,X]=multivariate_kernel_density(X0,N,kernel,id,findex,xindex,F,X)
%
%   [F,X]=multivariate_kernel_density(X0,N,kernel,id,findex,xindex,F,X,Hi_dHi)
%
% Inputs
% -------
%
% - **X0** [cell|matrix]: when it is a cell, the first element in the cell
% contains the observations and the second element contains the a d x 2
% matrix in which the first column represents lower bounds and the second
% columns represents upper bounds for the d variables of interest. When it
% is a matrix, X0 is simply a d x n matrix of observations, where d is the
% number of variables and n the number of observations. In that case,
% bounds default to the min and max in X0.
%
% - **N** [empty|integer|vector|{250}]: number of grid points in each
% dimension.
%
% - **kernel** [empty|char|{'normal'}]: must be one of the following:
%   'cosine','epanechnikov','normal','triangular','triweight','uniform'
%
% - **id** [empty|scalar]: indicates the current evaluation stage of the
% multivariate algorithm
%
% - **findex** [cell|empty]: coordinates to be filled in F below
%
% - **xindex** [empty|vector]: d-uple combination of grid points in each
% direction
%
% - **F** [empty|matrix]: kernel density
%
% - **X** [empty|cell]: each cell contains a 1 x N vector of coordinates
% for the observations in each dimension
%
% - **Hi_dHi** [cell]: placeholder for the inverse of the bandwith matrix
% and the determinant of the same inverse.
%
% Outputs
% --------
%
% - **F** [empty|matrix]: kernel density
%
% - **X** [empty|cell]: each cell contains a 1 x N vector of coordinates
% for the observations in each dimension
%
% More About
% ------------
%
% - In the two-dimensional case (d=2), the output is readily plotted using
% either mesh, surf or any suitable 3-D plotting function as:
%   - surf(X{1},X{2},F)
%   - mesh(X{1},X{2},F)
%   
% Examples
% ---------
%
% See also: kernel_density

% Sources: 
% - Wikipedia : 3.6 Multivariate Kernel Density Estimation
% - David W. Scott, Stephan R. Sain (): "Multi-dimensional Density
% Estimation", http://www.stat.rice.edu/~scottdw/ss.nh.pdf

if ~iscell(X0)
    lb=min(X0,[],2);
    ub=max(X0,[],2);
    X0={X0
        [lb,ub]};
end

[d,n]=size(X0{1});
if isscalar(N) && d>1
    N=N*ones(1,d);
end
if numel(N)~=d
    error('number of elements in N inconsistent with the number of variables')
end
N=N(:).';
matlab_plot_ready=true;
kernel_types={'cosine','epanechnikov','normal','triangular','triweight',...
    'uniform'};
if nargin<9
    Hi_dHi=[];
    if nargin<8
        X=[];
        if nargin<7
            F=[];
            if nargin<6
                xindex=[];
                if nargin<5
                    findex=[];
                    if nargin<4
                        id=[];
                        if nargin<3
                            kernel=[];
                        end
                    end
                end
            end
        end
    end
end

populate()

xlow=X0{2}(id,1);
xhigh=X0{2}(id,2);

X{id}=linspace(xlow,xhigh,N(id));

for ii=1:N(id)
    findex{id}=ii;
    xindex(id)=X{id}(ii);
    if id==d
        if matlab_plot_ready
            % somehow I need to do this so that things are not confused for
            % the surf and mesh functions later on
            flipped_findex=fliplr(findex);
            F(flipped_findex{:})=mvkd();
        else
            F(findex{:})=mvkd();
        end
    else
        [F,X]=multivariate_kernel_density(X0,N,kernel,id+1,findex,xindex,...
            F,X,Hi_dHi);
    end
end

    function populate()
        if isempty(Hi_dHi)
            %     h0=1;
            %     H=h0*eye(d);
            %     H=n^(-1/(d+4))*cov(X0.')^-.5;
            H=n^(-1/(d+4))*cov(X0{1}.')^.5;
            dHi=1/det(H);
            Hi=inv(H);
            Hi_dHi={Hi,dHi};
        end
        if isempty(X),X=cell(1,d); end
        if isempty(F)
            if d==1
                F=zeros(1,N);
            else
                if matlab_plot_ready
                    F=zeros(fliplr(N));
                else
                    F=zeros(N);
                end
            end
        end
        if isempty(xindex),xindex=zeros(d,1); end
        if isempty(findex),findex=cell(1,d); end
        if isempty(id),id=1; end
        if isempty(kernel),kernel='normal'; end
        if ~ismember(kernel,kernel_types)
            disp(kernel_types)
            error('kernel should be one of the above')
        end
    end
    function f=mvkd()
        xx_h=xindex(:,ones(1,n))-X0{1};
%         xx_h=bsxfun(@minus,xindex(:),X0{1});
        xx_h=Hi_dHi{1}*xx_h;
        %(1/n)*(1/det(H))*sum()
        f=1/n*Hi_dHi{2}*sum(kappa(xx_h));
        function k=kappa(u)
            nstar=cleanup();
            switch kernel
                case 'normal'
                    k=1/((2*pi)^(d/2))*exp(-0.5*sum(u.*u,1));
                case 'epanechnikov'
                    k=(3/4)^d*prod((1-u).^2,1);
                case 'triangular'
                    k=prod(1-abs(u));
                case 'triweight'
                    k=(35/32)^d*prod((1-u.^2).^3,1);
                case 'uniform'
                    k=(1/2)^d*ones(1,nstar);
                case 'cosine'
                    k=(pi/4)^d*prod(cos((pi/2)*u),1);
                otherwise
            end
            function nstar=cleanup()
                if ~strcmp(kernel,'normal')
                    bad=any(abs(u)>1,1);
                    u(:,bad)=[];
                end
                nstar=size(u,2);
            end
            
        end
    end
end
