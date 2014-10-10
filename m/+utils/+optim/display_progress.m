function display_progress(restart,iterations,fmin_global,fmin_iter,dispersion,funcCount,optimizer)
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
% Description
% ------------
%
% Examples
% ---------
%
% See also: 

if restart>1
	fprintf(1,'restart # %3.0f   iter: %6.0f   fmin(global) %8.4f    fmin(iter) %8.4f    stdev %8.4f    f-Count  %8.0f   routine %s\n',...
	        restart,iterations,fmin_global,fmin_iter,dispersion,funcCount,optimizer);
else
	fprintf(1,'iter: %6.0f   fmin(global) %8.4f    fmin(iter) %8.4f    stdev %8.4f    f-Count  %8.0f   routine %s\n',...
	        iterations,fmin_global,fmin_iter,dispersion,funcCount,optimizer);
end
		
% 1- restart
% 2- iterations
% 3- fmin(global)
% 4- fmin(iterations)
% 5- dispersion
% 6- funcCount
% 7- optimizer
