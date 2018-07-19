function RepsRun=bootstrap(self,boot)
% Bootstrap the parameter values resampling :math:`varepsilon_t` terms from the data sample
%
% ::
%
%    RepsRun = bootstrap(rfvar, boot)
%
% Args:
%    rfvar (rfvar object): rfvar object
%    boot (integer): number of bootstrap simulations (default: 1000)
%
% Returns:
%    :
%
%    - **RepsRun** : bootstrapped values of the estimated parameters
%
%       - dimension 1: parameter levels
%       - dimension 2: different runs of bootstrap samples
%
% Example:
%    For example, one can compute variance of the estimate of the parameters by computing::
%
%       RepsRun = bootstrap(rfvar, 10000);
%       var(RepsRun, [], 2)
%

if nargin < 2

    boot=[];

end

if isempty(boot),boot=1000; end

if ~isempty(self.estim_.prior)

    warning('No bootstrapping for models estimated using Bayesian techniques')

    if (self.is_switching||self.optimize)

        error('no bootstrapping for models estimated using bayesian techniques or optimzed numerically')

    else

        warning('switching to analytical sampling algorithms for constant-parameter BVARs')

        RepsRun=self.estim_.sampler(boot);

    end

end

sol=solve(self);

B=sol.B;

Resids=self.estim_.Y(:,:)-B*self.estim_.X(:,:);

% nv=self.nvars;
%
% nt=self.estim_.T;
%
% ng=self.ng;

Y=self.estim_.Y;

X=self.estim_.X;

T=self.estim_.T;

% K=kdata.K;

ndet=self.nx;

nvars=size(Y,1);

population=1:T;

fprintf('Starting bootstrapping %s\n', datestr(now));

for ii =1:boot

    SPoints = randsample(population,T,true);

    if ii==1

        RepsRun = cvarreps(Resids(:,SPoints));

        RepsRun=RepsRun(:,ones(boot,1));

    else

        RepsRun(:,ii) = cvarreps(Resids(:,SPoints));

    end

end

fprintf('bootstrapping reps finished, %s\n', datestr(now));

    function RepsRun = cvarreps(Resids)

        Yi=Y;

        Xi=X;

        for jcol=1:T

            Yi(:,jcol)=B*Xi(:,jcol)+Resids(:,jcol);

            if jcol<T

                Xi(ndet+1:end,jcol+1)=[Yi(:,jcol)
                    Xi(ndet+1:end-nvars,jcol)];

            end

        end

        kdatai=update_base_data(self,Xi,Yi);

        RepsRun = vartools.ols(kdatai);

        RepsRun=[RepsRun.B(:);vech(RepsRun.Sigma)];

        function k=update_base_data(k,X,Y)

            k.estim_.X=X;

            k.estim_.Y=Y;

        end

    end

end
