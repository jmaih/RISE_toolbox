function y1=one_step_engine(T,y0,ss,xloc,sig,shocks,order)

% H1 line
%
% ::
%
%
% Args:
%
% Returns:
%    :
%
% Note:
%
%    - It is expected that the solution of the first order also contains the
%    exogenous growth rate, which is nonzero for nonstationary models. That
%    component is the imaginary part of the column numel(xloc)+1.
%
% Example:
%
%    See also:

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

xloc0=xloc;

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

nlags=size(y0.y,2);

is_var=nlags>1;

pruned=isfield(y0,'y_lin');

% initialize output
%-------------------
y1=y0;

% initial conditions
%--------
if pruned
    
    y0=y0.y_lin;
    
    if isempty(y0)
        
        y0=y1.y;
        
    end
    
else
    
    y0=y0.y;
    
end

[n,no]=size(y0);

% demean
%--------
y0=y0-ss(:,ones(1,no));

% resize
%--------
if no<order
    
    if no==1
        
        y0=[y0,zeros(n,order-1)];
    
    else
        
        error('number of columns of y0 inconsistent with order of approximation')
    
    end
    
elseif (order<no||is_var) && order~=1
    
    error('order appears to be inconsistent with # columns of y0')
    
end

% extract the exogenous constant
%--------------------------------
try
    
    const_pos=numel(xloc)+1;
    
    Tsig_growth=T{1}(:,const_pos);
    
catch
    
    xloc=xloc0;
    
    const_pos=numel(xloc)+1;
    
    Tsig_growth=T{1}(:,const_pos);
    
end

growth=imag(Tsig_growth);

T{1}(:,const_pos)=real(Tsig_growth);

y1.y=0;

y01=zeros(n,order);

zkz=1;

if ~pruned
    
    z=stateify(1);
    
end

% add the various orders of approximation
%----------------------------------------
ifact=1./cumprod(1:order);

for iy=1:order
    
    for io=1*pruned+(1-pruned)*iy:iy
        
        zkz=mykron(io,zkz);
        
        y01(:,iy)=y01(:,iy)+ifact(io)*T{io}*zkz;
        
    end
    % cumulate main output
    %---------------------
    y1.y=y1.y+y01(:,iy);
    % now add the mean and store
    %----------------------------
    y01(:,iy)=y01(:,iy)+ss;
    
end

y1.y=y1.y+ss;
% add growth
y1.y=y1.y+growth;

if pruned
    
    y1.y_lin=y01;
    
    if any(growth)
        
        y1.y_lin(:,1)=y1.y_lin(:,1)+growth;
        
    end
    
end

    function zkz=mykron(io,zkz)
        
        if pruned
            
            zkz=0;
            
            [C,nrows]=find_combinations();
            
            for irow=1:nrows
                
                thisrow=C(irow,:);
                
                zi=stateify(thisrow(1));
                
                for icol=2:io
                    
                    zi=kron(zi,stateify(thisrow(icol)));
                    
                end
                
                zkz=zkz+zi;
                
            end
            
        else
            
            zkz=kron(zkz,z);
            
        end
        
        function [C,nrows]=find_combinations()
            % find all combinations (with repetition) of io elements such
            % that the sum is iy
            %--------------------------------------------------------------
            C=utils.gridfuncs.mygrid(iy*ones(1,io));
            
            test=sum(C,2)==iy;
            
            C=C(test,:);
            
            nrows=size(C,1);
            
        end
        
    end

    function z=stateify(id)
        
        if is_var
            
            x=y0(xloc,end:-1:1);% flip around if we have many lags
            
            x=x(:);
        
        else
            
            x=y0(xloc,id);
            
        end
        % nullify sig and shocks when id>1
        %---------------------------------
        zero_coef=1-(id>1);
        % form the state vector
        %-----------------------
        z=[x
            zero_coef*sig
            zero_coef*shocks(:)];
        
    end

end