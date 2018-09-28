function J=evaluate_jacobian_numerically(funcs,y,x,ss,param,def,s0,s1,param_id)
% INTERNAL FUNCTION: numerical evaluation of the jacobian of the objective function
%
% ::
%
%   J=evaluate_jacobian_numerically(funcs,y,x,ss,param,def,s0,s1)
%
% Args:
%
%    - **funcs** [fhandle|cell array]: function or functions to be
%      differentiated
%    - **y** [vector]: values of endogenous variables
%    - **x** [vector]: values of exogenous variables
%    - **ss** [vector]: steady state
%    - **param** [vector]: parameter vector
%    - **def** [vector]: values of definitions
%    - **s0** [scalar]: state today
%    - **s1** [scalar]: state tomorrow
%
% Returns:
%    :
%
%    - **J** [matrix]: Numerical Jacobian of **funcs** at [y,x,sparam]
%
% Note:
%    It is assumed that the inputs are y,x,ss,param,sparam,def,s0,s1 as
%    ordered in parser.input_list()
%

if nargin<9
    
    param_id=[];
    
end

ny=size(y,1);

nx=size(x,1);

np=numel(param_id);

yxp=[y;x;param(param_id)];%

if ~iscell(funcs)
    
    funcs={funcs};
    
end

[nrows,ncols]=size(funcs);

if ncols>1
    
    error('the objective cannot have many columns since outputs are to be concatenated vertically')
    
end

J=cell(nrows,1);

for irow=1:nrows
    
    J{irow}=utils.numdiff.jacobian(@newobjective,yxp);
    
end

if nrows==1
    
    J=J{1};
    
else
    
    J=cell2mat(J);
    
end

    function f=newobjective(yx)
        
        yy=yx(1:ny,:);
        
        xx=yx(ny+(1:nx),:);
        
        p=param;
        
        if np
            
            p(param_id)=yx(ny+nx+1:end);
            
        end
        
        f=funcs{irow}(yy,xx,ss,p,def,s0,s1);
        
    end

end