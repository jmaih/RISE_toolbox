function opts=update_options(TolFun,maxiter,debug)
% INTERNAL FUNCTION
%

opts=struct();

opts.fix_point_TolFun=TolFun;

opts.fix_point_maxiter=maxiter;

opts.fix_point_verbose=debug;

opts.fix_point_valid_func=@utils.error.validComplex;

end