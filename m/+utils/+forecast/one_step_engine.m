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