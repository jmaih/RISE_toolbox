function [sm]=pull_first_order_partitions(dv,posv)
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

sm=struct();

h=size(dv,2);

sm.dbf_plus=cell(h);

sm.ds_0=cell(h);

sm.dp_0=cell(h);

sm.db_0=cell(h);

sm.df_0=cell(h);

sm.dpb_minus=cell(h);

sm.de_0=cell(h);

for r0=1:h
    
    for r1=1:h
        
        sm.dbf_plus{r0,r1}=dv{r0,r1}(:,posv.bf_plus);
        
        sm.ds_0{r0,r1}=dv{r0,r1}(:,posv.s_0);
        
        sm.dp_0{r0,r1}=dv{r0,r1}(:,posv.p_0);
        
        sm.db_0{r0,r1}=dv{r0,r1}(:,posv.b_0);
        
        sm.df_0{r0,r1}=dv{r0,r1}(:,posv.f_0);
        
        sm.dpb_minus{r0,r1}=dv{r0,r1}(:,posv.pb_minus);
        
        sm.de_0{r0,r1}=dv{r0,r1}(:,posv.e_0);
    
    end
    
end

end
