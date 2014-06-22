function [T,itercode,retcode]=fix_point_iterator(iterate_func,T0,options,varargin) 
% this function solves for a fix point. Inputs are:
% 1-iterate_func: a function handle which returns 2 elements: The improved
% value of T and F0, the value of the function whose discrepancy is to be
% minmized 
% 2- T0 is the initial guess for the solution
% 3- options is a structure with fields 'fix_point_explosion_limit' (maximum
% divergence of F0), 'fix_point_TolFun', 'fix_point_maxiter', 'fix_point_verbose' 

default_solve=struct('fix_point_explosion_limit',1e+6,...
    'fix_point_TolFun',sqrt(eps),...
    'fix_point_maxiter',1000,...
    'fix_point_verbose',false);
if nargin==0
    if nargout>1
        error([mfilename,':: when the object is emtpy, nargout must be at most 1'])
    end
    T=default_solve;
    return
end

options=utils.miscellaneous.setfield(default_solve,options);
% solving for T
F0=0.1*options.fix_point_explosion_limit;
iter=0;
conv_T=inf;
conv_F=F0;
T00=T0;
fix_point_verbose=options.fix_point_verbose;
while max(conv_T,conv_F)>options.fix_point_TolFun && iter<options.fix_point_maxiter && ...
        conv_F<options.fix_point_explosion_limit && utils.error.valid(F0)
    iter=iter+1;
    [T,F0] = iterate_func(T0,varargin{:});
    conv_T=max(abs(T(:)-T0(:)));
    conv_F=max(abs(F0(:)));
    if fix_point_verbose
        fprintf(1,'%8.0f %12.4f %12.4f\n',iter,full(conv_T),conv_F); 
    end
    T0=T;
    if ~utils.error.valid(T)
        break
    end
end
retcode=0;
itercode=iter;
if iter>=options.fix_point_maxiter
    retcode=21;
elseif ~utils.error.valid(F0)||~utils.error.valid(T)
    retcode=22;
elseif conv_F>=options.fix_point_explosion_limit
    retcode=23;
end

if retcode==0 && ~utils.error.valid(T)
    assignin('base','T0',T00)
    assignin('base','iterate_func',iterate_func)
    assignin('base','options_',options)
    assignin('base','varargin_',varargin)
    disp(iter)
    keyboard
end


