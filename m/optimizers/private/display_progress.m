function display_progress(restart,iter,fmin_global,fmin_iter,dispersion,fcount,optimizer)

if restart>1
	fprintf(1,'restart # %3.0f   iter: %6.0f   fmin(global) %8.4f    fmin(iter) %8.4f    dispersion %8.4f    fcount  %8.0f   routine %s\n',...
	        restart,iter,fmin_global,fmin_iter,dispersion,fcount,optimizer);
else
	fprintf(1,'iter: %6.0f   fmin(global) %8.4f    fmin(iter) %8.4f    dispersion %8.4f    fcount  %8.0f   routine %s\n',...
	        iter,fmin_global,fmin_iter,dispersion,fcount,optimizer);
end
		
% 1- restart
% 2- iter
% 3- fmin(global)
% 4- fmin(iter)
% 5- dispersion
% 6- fcount
% 7- optimizer
