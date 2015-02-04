function [y1,iter,retcode]=one_step(T,y0,ss,xloc,sig,shocks,order,compl,cond_shocks_id)

% H1 line
%
% Syntax
% -------
% ::
%
% Inputs
% -------
%
% Outputs
% --------
%
% More About
% ------------
%
% Examples
% ---------
%
% See also:

narginchk(7,9)

n=nargin;
if n==7||isempty(compl)
    compl=@(x)1;
end
if n==7
    cond_shocks_id=[];
end
[nx,kplus1]=size(shocks);
if isempty(cond_shocks_id)
    cond_shocks_id=1:nx;
end

myzero=-sqrt(eps);
badvector=@(x)any(compl(x)<myzero);
fsolve_options=struct('Display','none','TolFun',1e-12,'TolX',1e-12);

y1=one_step_engine(T,y0,ss,xloc,sig,shocks,order);

iter=0;
retcode=0;
if badvector(y1.y)
    nrows=size(compl(y1.y),1);
    if islogical(cond_shocks_id)
        cond_shocks_id=find(cond_shocks_id);
    end
    n_cond_shocks=numel(cond_shocks_id);
    if n_cond_shocks~=nrows
        error('Singularity in the problem of matching restrictions')
    end
    y1k=y0(:,ones(1,kplus1+1));
    done=false;
    while ~done
        iter=iter+1;
        % find the shocks that make the violations go away and the
        % constraints bind exactly
        %------------------------------------------------------------------
        shocks0=zeros(nx,kplus1);
        shocks0(cond_shocks_id,1:iter)=nan;
        e_id=isnan(shocks0);
        ee0=zeros(n_cond_shocks*iter,1);
        [ee,fval,exitflag]=fsolve(@multi_complementarity,ee0,fsolve_options);
        exitflag=utils.optim.exitflag(exitflag,ee,max(abs(fval(:))));
        if exitflag~=1
            retcode=701;
            return
        end
        shocks0(e_id)=ee;
        % given those shocks, check that all future steps do not violate
        % the constraints
        %------------------------------------------------------------------
        test_passed=true;
        for icol=2:kplus1+1
            shocks_i=[shocks0(:,icol-1:end),zeros(nx,icol-2)];
            y1k(icol)=one_step_engine(T,y1k(icol-1),ss,xloc,sig,shocks_i,order);
            test_passed=~badvector(y1k(icol).y);
            if ~test_passed
                break
            end
            % here we could/should stop whenever there is a lift-up, since
            % there could be a chance that the step after the last one
            % still violates the constraints. But then this is true only if
            % the model is linear. If the model is nonlinear, it is
            % certainly a good idea to go to the end.
        end
        if iter==kplus1
            % the last guys is also constrained. Do one more step to see if
            % there is a lift up. If no lift up we are not out of the bush
            %--------------------------------------------------------------
            y1_last_uncond=one_step_engine(T,y1k(end),ss,xloc,sig,zeros(size(shocks0)),order);
            if ~all(compl(y1_last_uncond.y)>myzero)
                retcode=701;
                return
            end
        end
        % if the test is passed we are done
        %----------------------------------
        if test_passed
            y1=y1k(2);
            done=true;
        end
    end
end

    function viol=multi_complementarity(ee)
        shocks0(e_id)=ee;
        viol=zeros(nrows,iter);
        for icol_=2:iter+1
            shocks_i=[shocks0(:,icol_-1:end),zeros(nx,icol_-2)];
            y1k(icol_)=one_step_engine(T,y1k(icol_-1),ss,xloc,sig,shocks_i,order);
            viol(:,icol_-1)=compl(y1k(icol_).y);
        end
        viol=viol(:);
    end

end
%--------------------------------------------------------------------------
function y1=one_step_engine(T,y0,ss,xloc,sig,shocks,order)

% models with constants will have to be re-written as deviations from
% steady state. One could also have a zero steady state and instead put the
% constant as column in the shocks
if nargin<7
    order=numel(T);
end
if ~iscell(T)
    error('first input must be a cell')
end

if ~isstruct(y0)
    error('y0 should be a structure')
end

prunned=isfield(y0,'y_lin');
if isempty(xloc)
    % all variables are state variables
    xloc=1:size(y0.y,1);
end
if isempty(sig)
    % we assume by default that we do not have a model with sig. There will
    % be a crash if the matrices' sizes do not correspond to the state
    % vector.
    sig=0;
end

% build z
%--------
z=buildz(y0.y);
zkz=z;

% initialize y1.y at the steady state
%--------------------------------------
y1=y0;
y1.y=ss;

% add the various orders of approximation
%----------------------------------------
for io=1:order
    % account for VARs with many lags
    %--------------------------------
    y1.y=y1.y+1/factorial(io)*T{io}*zkz;
    if io==1 && prunned
        % vars will never go in here!
        if isempty(y0.y_lin)
            y1.y_lin=y1.y;
        else
            z_pruned=buildz(y0.y_lin);
            zkz=z_pruned;
            y1.y_lin=ss+1/factorial(io)*T{io}*zkz;
        end
    end
    if io<order
        if prunned && ~isempty(y0.y_lin)
            zkz=kron(zkz,z_pruned);
        else
            zkz=kron(zkz,z);
        end
    end
end

    function z=buildz(y00)
        % deviations from the steady state: use bsxfun in case we are in
        % presence of a VAR
        %---------------------------------------------------------------
        x=bsxfun(@minus,y00,ss);
        x=x(xloc,end:-1:1); % swap if we have many lags
        z=[
            x(:) % vectorize
            sig
            shocks(:)
            ];
        % take the sparse form to accelerate calculations in case of zero
        % shocks. This will also make the kroneckers sparse
        %----------------------------------------------------------------
        z=sparse(z);
    end
end