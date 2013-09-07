close all
figure
axis off
hold on
fill([-1 12 12 -1 -1],[-1 -1 12 12 -1],'w')
plot([-1 12 12 -1 -1],[-1 -1 12 12 -1],'k','Linewidth',3)
text(3,10,'Greek Symbols','color','r','FontSize',18)
x = 0; x1 = x+.7;
y = 4; y1 = y+.7;
z = 8; z1 = z+.7;
text(x,8,'\alpha :'); text(x1,8,'\\alpha');
text(x,7,'\beta :'); text(x1,7,'\\beta');
text(x,6,'\gamma :'); text(x1,6,'\\gamma')
text(x,5,'\delta :'); text(x1,5,'\\delta')
text(x,4,'\epsilon :'); text(x1,4,'\\epsilon')
text(x,3,'\kappa :'); text(x1,3,'\\kappa')
text(x,2,'\lambda :'); text(x1,2,'\\lambda')
text(x,1,'\mu :'); text(x1,1,'\\mu')
text(x,0,'\nu :'); text(x1,0,'\\nu')
text(y,8,'\omega :'); text(y1,8,'\\omega')
text(y,7,'\phi :'); text(y1,7,'\\phi')
text(y,6,'\pi :'); text(y1,6,'\\pi')
text(y,5,'\chi :'); text(y1,5,'\\chi')
text(y,4,'\psi :'); text(y1,4,'\\psi')
text(y,3,'\rho :'); text(y1,3,'\\rho')
text(y,2,'\sigma :'); text(y1,2,'\\sigma')
text(y,1,'\tau :'); text(y1,1,'\\tau')
text(y,0,'\upsilon :'); text(y1,0,'\\upsilon')
text(z,8,'\Sigma :'); text(z1,8,'\\Sigma')
text(z,7,'\Pi :'); text(z1,7,'\\Pi')
text(z,6,'\Lambda :'); text(z1,6,'\\Lambda')
text(z,5,'\Omega :'); text(z1,5,'\\Omega')
text(z,4,'\Gamma :'); text(z1,4,'\\Gamma')
shg
%%
greeks = ...
{'ALPHA'      'A'         '\alpha'
 'BETA'       'B'         '\beta'
 'GAMMA'      '\Gamma'    '\gamma'
 'DELTA'      '\Delta'    '\delta'
 'EPSILON'    'E'         {'\epsilon','\varepsilon'}
 'ZETA'       'Z'         '\zeta'
 'ETA'        'H'         '\eta'
 'THETA'      '\Theta'    {'\theta','\vartheta'}
 'IOTA'       'I'         '\iota'
 'KAPPA'      'K'         '\kappa'
 'LAMBDA'     '\Lambda'   '\lambda'
 'MU'         'M'         '\mu'
 'NU'         'N'         '\nu'
 'XI'         '\Xi'       '\xi'
 'OMICRON'    'O'         'o'
 'PI'         '\Pi'       {'\pi','\varpi'}
 'RHO'        'P'         {'\rho','\varrho'}
 'SIGMA'      '\Sigma'    {'\sigma','\varsigma'}
 'TAU'        'T'         '\tau'
 'UPSILON'    '\Upsilon'  '\upsilon'
 'PHI'        '\Phi'      {'\phi','\varphi'}
 'CHI'        'X'         '\chi'
 'PSI'        '\Psi'      '\psi'
 'OMEGA'      '\Omega'    '\omega'};
h = figure('units','pixels','pos',[300,100,620,620],'Color','w');
axes('units','pixels','pos',[10,10,600,600],'Xcol','w','Ycol','w',...
     'Xtick',[],'Ytick',[],'Xlim',[0 6],'Ylim',[0,4]);
% Loop by column and row

for r = 1:4
    for c = 1:6
        el = (r-1)*6 + c;
        % Title
        text(c-0.5,5-r,greeks{el,1},'Fonts',14,'FontN','FixedWidth',...
            'Hor','center','Ver','cap')
        % Color cap latter in grey or black
        if strcmp(greeks{el,2}(1),'\')
            clr = [0, 0, 0];
        else
            clr = [0.65, 0.65, 0.65];
        end
        % Cap letter
        text(c-0.5,4.87-r,['$\rm{' greeks{el,2} '}$'],'Fonts',40,...
            'Hor','center','Ver','cap','Interp','Latex','Color',clr)
        % Lowercase letter/s (if two variants) 
        if iscell(greeks{el,3})
            text(c-0.75,4.48-r,['$' greeks{el,3}{1} '$'],'Fonts',20,...
                'Hor','center','Interp','Latex')
            text(c-0.25,4.48-r,['$' greeks{el,3}{2} '$'],'Fonts',20,...
                'Hor','center','Interp','Latex')
            % Latex command
            text(c-0.5,4.3-r,['\' greeks{el,3}{1}],'Fonts',12,'FontN','FixedWidth',...
            'Hor','center','Ver','base')
        else
            text(c-0.5,4.48-r,['$' greeks{el,3} '$'],'Fonts',20,...
                'Hor','center','Interp','Latex')
            text(c-0.5,4.3-r,['\' greeks{el,3}],'Fonts',12,'FontN','FixedWidth',...
            'Hor','center','Ver','base')
        end
    end
end