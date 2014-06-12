%% housekeeping
clear classes
clc
%%
r=rise_svar.template();
r.model='svar';
r.constant=true;
r.nlags=3;
r.endogenous={
    'robs'   ,'interest rates'
    'dy'     ,'GDP growth'
    'labobs' ,'hours worked'
    'pinfobs','inflation'
    'dw'     ,'wages growth'
};
r.exogenous={
    'mp','Monetary Policy'
    'ad','Aggregate demand'
    'as','Aggregate supply'
    'ls','Labor supply'
    'junk','not important'
    };
r.short_run_restrictions={
    'dy@mp',0
    'dy@ad',0
    'dy@as',0
    'dy@ls',0
    };
r.long_run_restrictions={
    };
r.sign_restrictions={
    'robs@mp','+'
    'robs@ad','+'
    'robs@as','-'
    'robs@ls','-'
    'dy@mp','-'
    'dy@ad','+'
    'dy@as','+'
    'dy@ls','+'
    'labobs@ls','+'
    'pinfobs@mp','-'
    'pinfobs@ad','+'
    'pinfobs@as','-'
    'pinfobs@ls','-'
    'dw@ls','-'
};
%% Bring in the data
[num,txt,raw] = xlsread('Smets_Wouters_data.xlsx');

date_start='1947q3';
robs_id=2; dc_ic=3; dinve_id=4; dy_id=5; labobs_id=6; pinfobs_id=7; dw_id=8;
data=ts(date_start,...
    num(:,[robs_id,dy_id,labobs_id,pinfobs_id,dw_id]),...
    {'robs','dy','labobs','pinfobs','dw'});
data1=pages2struct(data);
%% construct the rise_svar object
rv=rise_svar(r,'data',data);
%%
rv=estimate_reduced_form(rv);
rv=find_structural_form(rv);
disp('========== A0 matrix ==============')
disp(rv.sf_A0{1})
disp('========== Impact matrix ==============')
disp([{[]},rv.exo_names;rv.endo_names',num2cell(inv(rv.sf_A0{1}))])
%% Plot some impulse responses
myirfs=irf(rv);
out=fanchart(myirfs.mp,[0.3,0.5,0.9]);
plot_fanchart(out.dy,'b')
%%
clc
A=cell(0,2);
%%
A0=nan(3); A0(1,2)=0; A0(2,3)=0; A0(3,1)=0;
A=[A;{A0,[]}];
%%
A0=nan(3); A0(1:2,1)=0; A0(3,2)=0;
A=[A;{A0,[]}];
%%
A0=nan(3); A0(1:2,1)=0; A0(3,2:3)=0; 
A1=nan(3); A1(1:2,1)=0; A1(3,2)=0;
A=[A;{A0,A1}];
%%
A0=nan(5);
A0(2:5,1)=0;A0(3:5,2)=0;A0([1,2,5],3)=0;A0(5,4)=0;
A=[A;{A0,[]}];
%% 
rv=rise_svar.empty(0);
for irow=1:size(A,1)
    r=rise_svar.template();
    A0=A{irow,1};
    A1=A{irow,2};
    n=size(A0,1);
    r.endogenous=cellstr(strcat('endo_',num2str((1:n)')));
    r.exogenous=cellstr(strcat('exo_',num2str((1:n)')));
    for irow_a0=1:size(A0,1)
        for icol_a0=1:n
            if A0(irow_a0,icol_a0)==0
                v=[sprintf('endo_%0.0f',irow_a0),'@',sprintf('exo_%0.0f',icol_a0)];
                r.short_run_restrictions=[r.short_run_restrictions
                    {v,0}];
            end
        end
    end
    for irow_a1=1:size(A1,1)
        for icol_a1=1:n
            if A1(irow_a1,icol_a1)==0
                v=[sprintf('endo_%0.0f',irow_a1),'{-1}@',sprintf('exo_%0.0f',icol_a1)];
                r.short_run_restrictions=[r.short_run_restrictions
                    {v,0}];
            end
        end
    end
    rv(irow,1)=rise_svar(r);
    fprintf(1,'============model %0.0f =========\n',irow);
    disp(A0)
    disp(A1)
    disp(rv(irow))
end