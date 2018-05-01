function [y1,is_active_shock,retcode]=one_step(T,y0,ss,xloc,sig,shocks,...
order,compl,shock_structure)

% H1 line
%
% ::
%
%
% Args:
%
% Returns:
%    :
%
% Note:
%
% Example:
%
%    See also:

narginchk(7,9)

n=nargin;
if n==7
    compl=[];
end
if n==7
    shock_structure=[];
end
[~,kplus1]=size(shocks);

myzero=-sqrt(eps);
% badvector=@(x)any(compl(x)<myzero); % recreating this slows things down...
debug=true;
if debug
    fsolve_options=struct('Display','iter','TolFun',1e-12,'TolX',1e-12);
else
    fsolve_options=struct('Display','none','TolFun',1e-12,'TolX',1e-12);
end

y1=utils.forecast.one_step_engine(T,y0,ss,xloc,sig,shocks,order);

retcode=0;
is_active_shock=false(1,kplus1);
if ~isempty(compl) && any(compl(y1.y)<myzero) % <---badvector(y1.y)
    if ~islogical(shock_structure)
        error('shock_structure should be logical');
    end
    y1k=y0(:,ones(1,kplus1+1));
    % find the shocks that make the violations go away and the
    % constraints bind exactly
    %------------------------------------------------------------------
    shocks0=shocks;
    % let the contemporaneous shocks be
    shocks0(shock_structure)=nan;
    e_id=isnan(shocks0);
    ee0=zeros(sum(e_id(:)),1);
    [ee,fval,exitflag]=fsolve(@multi_complementarity,ee0,fsolve_options);
    exitflag=utils.optim.exitflag(exitflag,ee,max(abs(fval(:))));
    if exitflag~=1
        retcode=701;
        return
    end
    shocks0(e_id)=ee;
    y1=y1k(2);
    is_active_shock(2:end)=true;
end

    function viol=multi_complementarity(ee)
        shocks0(e_id)=ee;
        y1k(2)=utils.forecast.one_step_engine(T,y1k(1),ss,...
            xloc,sig,shocks0,order);
        viol=compl(y1k(2).y);
    end

end