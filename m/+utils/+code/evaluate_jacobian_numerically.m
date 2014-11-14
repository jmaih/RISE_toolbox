function J=evaluate_jacobian_numerically(funcs,y,x,ss,param,sparam,def,s0,s1,switch_params_leads_index)
% evaluate_jacobian_numerically - numerical evaluation of the jacobian of the objective function
%
% Syntax
% -------
% ::
%
%   J=evaluate_jacobian_numerically(funcs,y,x,ss,param,sparam,def,s0,s1)
%
% Inputs
% -------
%
% - **funcs** [fhandle|cell array]: function or functions to be
%   differentiated
%
% - **y** [vector]: values of endogenous variables
%
% - **x** [vector]: values of exogenous variables
%
% - **ss** [vector]: steady state 
%
% - **param** [vector]: parameter vector 
%
% - **sparam** [vector]: vector of parameters appearing with a lead 
%
% - **def** [vector]: values of definitions 
%
% - **s0** [scalar]: state today 
%
% - **s1** [scalar]: state tomorrow 
%
% - **switch_params_leads_index** [[]|scalar|vector]: indices of the
%   switching parameters appearing with a lead
%
% Outputs
% --------
%
% - **J** [matrix]: Numerical Jacobian of **funcs** at [y,x,sparam]
%
% More About
% ------------
%
% It is assumed that the inputs are y,x,ss,param,sparam,def,s0,s1 as
% ordered in parser.input_list()
%
% Examples
% ---------
%
% See also: 

ny=size(y,1);
nx=size(x,1);
yx=[y;x;sparam(switch_params_leads_index)];%
nyxsp=numel(yx);
sp=sparam(:,ones(1,nyxsp));
if ~iscell(funcs)
    funcs={funcs};
end
[nrows,ncols]=size(funcs);
if ncols>1
    error('the objective cannot have many columns since outputs are to be concatenated vertically')
end

J=cell(nrows,1);
for irow=1:nrows
    J{irow}=utils.numdiff.jacobian(@newobjective,yx);
end
if nrows==1
    J=J{1};
else
    J=cell2mat(J);
end

    function f=newobjective(yx)
        yy=yx(1:ny,:);
        xx=yx(ny+(1:nx),:);
        sp(switch_params_leads_index,:)=yx(ny+nx+1:end,:);
        f=funcs{irow}(yy,xx,ss,param,sp,def,s0,s1);
    end
end