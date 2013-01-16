function display_progress(restart,iterations,fmin_global,fmin_iter,dispersion,funcCount,optimizer)

if restart>1
	fprintf(1,'restart # %3.0f   iterations: %6.0f   fmin(global) %8.4f    fmin(iterations) %8.4f    dispersion %8.4f    funcCount  %8.0f   routine %s\n',...
	        restart,iterations,fmin_global,fmin_iter,dispersion,funcCount,optimizer);
else
	fprintf(1,'iterations: %6.0f   fmin(global) %8.4f    fmin(iterations) %8.4f    dispersion %8.4f    funcCount  %8.0f   routine %s\n',...
	        iterations,fmin_global,fmin_iter,dispersion,funcCount,optimizer);
end
		
% 1- restart
% 2- iterations
% 3- fmin(global)
% 4- fmin(iterations)
% 5- dispersion
% 6- funcCount
% 7- optimizer
