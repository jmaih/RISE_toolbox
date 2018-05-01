function [F,X]=multivariate_kernel_density(X0,N,kernel,id,findex,F,X,Hi_dHiOvern)
% multivariate_kernel_density -- multivariate kernel density estimation
%
% ::
%
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
%   [F,X]=multivariate_kernel_density(X0,N,kernel,id,findex,F)
%
%   [F,X]=multivariate_kernel_density(X0,N,kernel,id,findex,F,X)
%
%   [F,X]=multivariate_kernel_density(X0,N,kernel,id,findex,F,X,Hi_dHi)
%
% Args:
%
%    - **X0** [cell|matrix]: when it is a cell, the first element in the cell
%    contains the observations and the second element contains the a d x 2
%    matrix in which the first column represents lower bounds and the second
%    columns represents upper bounds for the d variables of interest. When it
%    is a matrix, X0 is simply a d x n matrix of observations, where d is the
%    number of variables and n the number of observations. In that case,
%    bounds default to the min and max in X0.
%
%    - **N** [empty|integer|vector|{250}]: number of grid points in each
%    dimension.
%
%    - **kernel** [empty|char|{'normal'}]: must be one of the following:
%      'cosine','epanechnikov','normal','triangular','triweight','uniform'
%
%    - **id** [empty|scalar]: indicates the current evaluation stage of the
%    multivariate algorithm
%
%    - **findex** [cell|empty]: coordinates to be filled in F below
%
%    - **F** [empty|matrix]: kernel density
%
%    - **X** [empty|cell]: each cell contains a 1 x N vector of coordinates
%    for the observations in each dimension
%
%    - **Hi_dHiOvern** [cell]: placeholder for the inverse of the bandwith
%    matrix and the determinant of the same inverse divided by n.
%
% Returns:
%    :
%
%    - **F** [empty|matrix]: kernel density
%
%    - **X** [empty|cell]: each cell contains a 1 x N vector of coordinates
%    for the observations in each dimension
%
% Note:
%
%    - In the two-dimensional case (d=2), the output is readily plotted using
%    either mesh, surf or any suitable 3-D plotting function as:
%      - surf(X{1},X{2},F)
%      - mesh(X{1},X{2},F)
%
% Example:
%
%    See also: kernel_density

% Sources:
% - Wikipedia : 3.6 Multivariate Kernel Density Estimation
% - David W. Scott, Stephan R. Sain (): "Multi-dimensional Density
% Estimation", http://www.stat.rice.edu/~scottdw/ss.nh.pdf
% - wikipedia: https://en.wikipedia.org/wiki/Kernel_(statistics)

lb=[];
ub=[];
if iscell(X0)
    lb=X0{2}(:,1);
    ub=X0{2}(:,2);
    X0=X0{1};
end

[d,n]=size(X0);
if isscalar(N) && d>1
    N=N*ones(1,d);
end
if numel(N)~=d
    error('number of elements in N inconsistent with the number of variables')
end
N=N(:).';
matlab_plot_ready=true;
non_truncated_kernels={'normal','gaussian','logistic','silverman'};
kernel_types=[non_truncated_kernels,{'cosine','epanechnikov','biweight',...
    'quartic','triangular','triweight','tricube','uniform'}];
if nargin<8
    Hi_dHiOvern=[];
    if nargin<7
        X=[];
        if nargin<6
            F=[];
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

if isempty(id)
    populate()
end

X00=X0;
x0_id=X0(id,:);
x_id=X{id};
for ii=1:N(id)
    findex{id}=ii;
    X00(id,:)=x_id(ii)-x0_id;
    if id==d
        if matlab_plot_ready
            % somehow I need to do this so that things are not confused for
            % the surf and mesh functions later on
            flipped_findex=fliplr(findex);
            F(flipped_findex{:})=mvkd____();
        else
            F(findex{:})=mvkd____();
        end
    else
        [F,X]=multivariate_kernel_density(X00,N,kernel,id+1,findex,F,X,...
            Hi_dHiOvern);
    end
end

    function populate()
        % bandwith
        %------------
        %     h0=1;
        %     H=h0*eye(d);
        %     H=n^(-1/(d+4))*cov(X0.')^-.5;
        H=n^(-1/(d+4))*cov(X0.')^.5;
        dHi_n=1/det(H)*1/n;
        Hi=inv(H);
        Hi_dHiOvern={Hi,dHi_n};
        % grid points
        %------------
        if isempty(lb)
            lb=min(X0,[],2);
            ub=max(X0,[],2);
        end
        X=cell(1,d);
        for ivar=1:d
            X{ivar}=linspace(lb(ivar),ub(ivar),N(ivar));
        end
        % density
        %--------
        if d==1
            F=zeros(1,N);
        else
            if matlab_plot_ready
                F=zeros(fliplr(N));
            else
                F=zeros(N);
            end
        end
        % indexing of density
        %--------------------
        findex=cell(1,d);
        % variable index
        %----------------
        id=1;
        % kernel density type
        %---------------------
        if isempty(kernel),kernel='normal'; end
        if ~ismember(kernel,kernel_types)
            disp(kernel_types)
            error('kernel should be one of the above')
        end
    end
    function f=mvkd____()
        xx_h=Hi_dHiOvern{1}*X00;
        %(1/n)*(1/det(H))*sum()
        f=kappa(xx_h);
        f=Hi_dHiOvern{2}*sum(f);
        function k=kappa(u)
            if ~any(strcmp(kernel,non_truncated_kernels))
                nstar=[];
                cleanup()
            end
            switch kernel
                case 'uniform'
                    k=(1/2)^d*ones(1,nstar);
                case 'triangular'
                    k=1-abs(u);
                    k=prod(k,1);
                case 'epanechnikov'
                    k=(3/4)*(1-u.^2);
                    k=prod(k,1);
%                 case 'spherical_epanechnikov'
%                     prop=1;
%                     k=sum(u.*u,1);
%                     k=k(k<=1);
%                     k=prop*(1-k);
                case {'quartic','biweight'}
                    k=(15/16)*(1-u.^2).^2;
                    k=prod(k,1);
                case 'triweight'
                    k=(35/32)*(1-u.^2).^3;
                    k=prod(k,1);
                case 'tricube'
                    k=(70/81)*(1-u.^3).^3;
                    k=prod(k,1);
                case {'normal','gaussian'}
                    k=1/((2*pi)^(0.5*d))*exp(-0.5*sum(u.*u,1));
                case 'cosine'
                    k=(0.25*pi)*cos(0.5*pi*u);
                    k=prod(k,1);
                case 'logistic'
                    k=1./(exp(u)+2+exp(-u));
                    k=prod(k,1);
                case 'silverman'
                    k=abs(u)/sqrt(2);
                    k=0.5*exp(-k)*sin(k+0.25*pi);
                    k=prod(k,1);
            end
            function cleanup()
                bad=any(abs(u)>1,1);
                u(:,bad)=[];
                nstar=size(u,2);
            end
        end
    end
end
