function [ss,newp,retcode]=ssfile(m,ss,p,d,id)
% [y,newp,retcode]=steady_state_4_Canonical_Const(obj,y,p,d,id)

retcode=0;
if nargin==1
    % list of endogenous variables to be calculated
    %----------------------------------------------
    ss={'aj','ak','ap','arr','aw','az','b','bnot','c','c1','dp','dp1',...
        'dp2','dp3','dw','dw1','h1','ik','k','lm','n','n1','q','qk','r',...
        'rk','rnot','uc','uc1','uh','uh1','w','w1','xp','xw','xw1','y'};
    % list of parameters to be computed during steady state calculation
    %-------------------------------------------------------------------
    newp={};
    
    return
    
end

newp=[];
    
r = p.PIBAR / p.BETA ;
rk = 1/p.BETA - (1-p.DK) ;
xp = p.XP_SS ;
xw = p.XW_SS ;
xw1 = p.XW_SS ;
lm = (1 - p.BETA1/p.BETA) / (1 - p.BETA1*p.RHOD/p.PIBAR) ;

QHTOC_tmp = p.JEI/(1-p.BETA);
QH1TOC1_tmp = p.JEI/(1-p.BETA1-lm*(1-p.RHOD)*p.M);
KTOY_tmp = p.ALPHA/(xp*rk);
BTOQH1_tmp = p.M*(1-p.RHOD)/(1-p.RHOD/p.PIBAR) ;
C1TOY_tmp = (1-p.ALPHA)*p.SIGMA/(1+(1/p.BETA-1)*BTOQH1_tmp*QH1TOC1_tmp)*(1/xp) ;
CTOY_tmp = (1-C1TOY_tmp-p.DK*KTOY_tmp) ;

n = ((1-p.SIGMA)*(1-p.ALPHA)/(xp*xw*CTOY_tmp))^(1/(1+p.ETA));
n1 = (p.SIGMA*(1-p.ALPHA)/(xp*xw1*C1TOY_tmp))^(1/(1+p.ETA));

y = KTOY_tmp^(p.ALPHA/(1-p.ALPHA))*(n^(1-p.SIGMA))*n1^p.SIGMA ;

c = CTOY_tmp*y;
c1 = C1TOY_tmp*y;
k = KTOY_tmp*y;
ik = p.DK*k;

w = xw*c*n^p.ETA;
w1 = xw1*c1*n1^p.ETA;
q = QHTOC_tmp*c + QH1TOC1_tmp*c1 ;
h = QHTOC_tmp*c/q ;
h1 = QH1TOC1_tmp*c1/q ;
b = BTOQH1_tmp*q*h1 ;
uc = 1/c;
uc1 = 1/c1;
uh = p.JEI/h;
uh1 = p.JEI/h1;
dp = p.PIBAR ;
dp1 = dp;
dp2 = dp;
dp3 = dp;
dw = p.PIBAR ;
dw1 = p.PIBAR ;
aj = 1 ;
arr =1;
az = 1;
ak=1;
ap=1;
aw=1;
% un = n^ETA ;
% un1 = n1^ETA ;
% aa = 1;
% af=1;
% am = 1;
% an=1;

qk = 1;
rnot = r;

% lev = b/(q*h1);
bnot = b;
  

ys = [aj ak ap arr aw az b bnot c c1 dp  dp1 dp2 dp3 dw dw1 h1 ik k lm n ...
    n1 q qk r rk rnot uc uc1 uh uh1 w w1 xp xw xw1 y ];

ss(id)=ys;


end


