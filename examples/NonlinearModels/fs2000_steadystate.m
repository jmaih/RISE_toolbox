% computes the steady state of fs2000 analyticaly
function [ys,retcode]=fs2000_steadystate(param_obj,flag)

retcode=0;
switch flag
    case 0
        ys={'m','P','c','e','W','R','k','d','n','l','gy_obs','gp_obs','y','dA'};
    case 1
        param_names={param_obj.name};
        params=vertcat(param_obj.startval); %#ok<NASGU>
        
        for index=1:numel(param_names)
            eval([param_names{index},'=params(index);'])
        end
        
        % N.B: psi is the name of a function in Matlab. For some reason I do not
        % understand, it creates problems when running this file although the loop
        % above creates a parameter psi and assigns it a value...
        
        dA = exp(gam);
        gst = 1/dA;
        m = mst;
        
        khst = ( (1-gst*bet*(1-del)) / (alp*gst^alp*bet) )^(1/(alp-1));
        xist = ( ((khst*gst)^alp - (1-gst*(1-del))*khst)/mst )^(-1);
        eval('nust = psi*mst^2/( (1-alp)*(1-psi)*bet*gst^alp*khst^alp );')
        n  = xist/(nust+xist);
        P  = xist + nust;
        k  = khst*n;
        
        eval('l  = psi*mst*n/( (1-psi)*(1-n) );')
        c  = mst/P;
        d  = l - mst + 1;
        y  = k^alp*n^(1-alp)*gst^alp;
        R  = mst/bet;
        W  = l/n;
        %   ist  = y-c;
        %   q  = 1 - d;
        
        e = 1;
        
        gp_obs = m/dA;
        gy_obs = dA;
        
        ys =[m P c e W R k d n l gy_obs gp_obs y dA]';
    otherwise
        error([mfilename,':: Unknown flag'])
end

