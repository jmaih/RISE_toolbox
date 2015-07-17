function [dbf_plus,ds_0,dp_0,db_0,df_0,dpb_minus,de_0]=pull_first_order_partitions(dv,posv)
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
% More About
% ------------
%
% Examples
% ---------
%
% See also: 

h=size(dv,2);
dbf_plus=cell(h);
ds_0=cell(h);
dp_0=cell(h);
db_0=cell(h);
df_0=cell(h);
dpb_minus=cell(h);
de_0=cell(h);
for r0=1:h
    for r1=1:h
        dbf_plus{r0,r1}=dv{r0,r1}(:,posv.bf_plus);
        ds_0{r0,r1}=dv{r0,r1}(:,posv.s_0);
        dp_0{r0,r1}=dv{r0,r1}(:,posv.p_0);
        db_0{r0,r1}=dv{r0,r1}(:,posv.b_0);
        df_0{r0,r1}=dv{r0,r1}(:,posv.f_0);
        dpb_minus{r0,r1}=dv{r0,r1}(:,posv.pb_minus);
        de_0{r0,r1}=dv{r0,r1}(:,posv.e_0);
    end
end
end
