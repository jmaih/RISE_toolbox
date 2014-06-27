% computes the steady state of fs2000 analyticaly
function [ys,param_obj,retcode,imposed]=fs2000_steadystate(param_obj,flag)

retcode=0;

imposed=false;
% setting this to false tells rise that this is just an initial guess and so, rise
% will check that it is actually the steady state and if it is not, rise will use
% it as an initial value
switch flag
    case 0
        ys={'m','P','c','e','W','R','k','d','n','l','gy_obs','gp_obs','y','dA'};
    case 1
        pp=struct();
        name_loc=strcmp('name',param_obj(1,:));
        val_loc=strcmp('startval',param_obj(1,:));
        par_names=param_obj{2,name_loc};
        par_mat=param_obj{2,val_loc};
        for ipar=1:numel(par_names)
            pp.(par_names{ipar})=par_mat(ipar,1);
        end
                
        dA = exp(pp.gam);
        gst = 1/dA;
        m = pp.mst;
        
        khst = ( (1-gst*pp.bet*(1-pp.del)) / (pp.alp*gst^pp.alp*pp.bet) )^(1/(pp.alp-1));
        xist = ( ((khst*gst)^pp.alp - (1-gst*(1-pp.del))*khst)/pp.mst )^(-1);
        nust = pp.psi*pp.mst^2/( (1-pp.alp)*(1-pp.psi)*pp.bet*gst^pp.alp*khst^pp.alp );
        n  = xist/(nust+xist);
        P  = xist + nust;
        k  = khst*n;
        
        l  = pp.psi*pp.mst*n/( (1-pp.psi)*(1-n) );
        c  = pp.mst/P;
        d  = l - pp.mst + 1;
        y  = k^pp.alp*n^(1-pp.alp)*gst^pp.alp;
        R  = pp.mst/pp.bet;
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

