function [ys,retcode]=tshocksnk_steadystate(param_obj,flag)

retcode=0;
switch flag
    case 0
        ys={'ln_pai','ln_a','ln_z','ln_theta','ln_r','ln_g','ln_c','ln_y','ln_x'};
    case 1
        param_names={param_obj.name};
        params=vertcat(param_obj.startval); %#ok<NASGU>
        
        % using greek letters that also functions in matlab seems to
        % consistently cause problems. One way to deal with the problem is
        % to initialize the variable....
        beta=[];
        for index=1:numel(param_names)
            eval([param_names{index},'=params(index);'])
        end
        eta=1/omega;
%         phi=eta*(theta_ss-1)/psi;

		pai=pai_ss;
		a=a_ss;
		z=z_ss;
		theta=theta_ss;
		r=pai*z/beta;
		g=z;
		c=(a*((theta-1)/theta))^(1/eta);
		y=c;
		x=((theta-1)/theta)^(1/eta);
        
        ys =log([pai,a,z,theta,r,g,c,y,x])';
    otherwise
        error([mfilename,':: Unknown flag'])
end


