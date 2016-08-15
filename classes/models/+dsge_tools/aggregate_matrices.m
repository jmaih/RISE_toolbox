function [sm,adjusted,accelSupport]=aggregate_matrices(sm,siz,adjusted,accelerate)

accelSupport=struct();

% aggregate A0 and A_
%--------------------
sm.d0=cell(1,siz.h);

for r0=1:siz.h
    
    for r1=2:siz.h
        
        sm.ds_0{r0,1}=sm.ds_0{r0,1}+sm.ds_0{r0,r1};
        
        sm.dp_0{r0,1}=sm.dp_0{r0,1}+sm.dp_0{r0,r1};
        
        sm.db_0{r0,1}=sm.db_0{r0,1}+sm.db_0{r0,r1};
        
        sm.df_0{r0,1}=sm.df_0{r0,1}+sm.df_0{r0,r1};
        
        sm.dpb_minus{r0,1}=sm.dpb_minus{r0,1}+sm.dpb_minus{r0,r1};
        
        sm.de_0{r0,1}=sm.de_0{r0,1}+sm.de_0{r0,r1};
        
    end
    
    sm.d0{r0}=[sm.ds_0{r0,1},sm.dp_0{r0,1},sm.db_0{r0,1},sm.df_0{r0,1}];
    
    % eliminate static variables for speed
    %-------------------------------------
    if accelerate
        
        if r0==1
            
            accelSupport.Abar_minus_s=cell(1,siz.h);
            
            accelSupport.R_s_s=cell(1,siz.h);
            
            accelSupport.R_s_ns=cell(1,siz.h);
            
            accelSupport.Abar_plus_s=cell(siz.h);
            
            adjusted.nd=adjusted.nd-siz.ns;
            
            adjusted.bf_cols=adjusted.bf_cols-siz.ns;
            
            adjusted.pb_cols=adjusted.pb_cols-siz.ns;
            
            adjusted.siz.ns=0;
            
            adjusted.siz.nd=adjusted.siz.nd-siz.ns;
            
            adjusted.siz.nT=adjusted.siz.nT-siz.ns;
            
        end
        
        [Q0,sm.d0{r0}]=qr(sm.d0{r0});
        
        sm.dpb_minus{r0,1}=Q0'*sm.dpb_minus{r0,1};
        
        accelSupport.Abar_minus_s{r0}=sm.dpb_minus{r0,1}(1:siz.ns,:);
        
        sm.dpb_minus{r0,1}=sm.dpb_minus{r0,1}(siz.ns+1:end,:);
        
        for r1=1:siz.h
            
            sm.dbf_plus{r0,r1}=Q0'*sm.dbf_plus{r0,r1};
            
            accelSupport.Abar_plus_s{r0,r1}=sm.dbf_plus{r0,r1}(1:siz.ns,:);
            
            sm.dbf_plus{r0,r1}=sm.dbf_plus{r0,r1}(siz.ns+1:end,:);
            
        end
        
        accelSupport.R_s_s{r0}=sm.d0{r0}(1:siz.ns,1:siz.ns);
        
        accelSupport.R_s_ns{r0}=sm.d0{r0}(1:siz.ns,siz.ns+1:end);
        
        sm.d0{r0}=sm.d0{r0}(siz.ns+1:end,siz.ns+1:end);
        
        sm.ds_0{r0,1}=sm.d0{r0}(:,1:adjusted.siz.ns);
        
        sm.dp_0{r0,1}=sm.d0{r0}(:,adjusted.siz.ns+(1:adjusted.siz.np));
        
        sm.db_0{r0,1}=sm.d0{r0}(:,adjusted.siz.ns+adjusted.siz.np+(1:adjusted.siz.nb));
        
        sm.df_0{r0,1}=sm.d0{r0}(:,adjusted.siz.ns+adjusted.siz.np+adjusted.siz.nb+(1:adjusted.siz.nf));
        
    end
    
end

sm.ds_0=sm.ds_0(:,1)';

sm.dp_0=sm.dp_0(:,1)';

sm.db_0=sm.db_0(:,1)';

sm.df_0=sm.df_0(:,1)';

sm.de_0=sm.de_0(:,1)';

end