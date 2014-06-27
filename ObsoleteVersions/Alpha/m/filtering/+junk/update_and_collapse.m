function  [att,Ptt,K]=update_and_collapse(a01,P01,iF,v,obs_id,endo_nbr,pp,h,...
    hstar,PAI01_tt,PAItt)
a01tt=zeros(endo_nbr,h,hstar);
P01tt=zeros(endo_nbr,endo_nbr,h,hstar);
att=zeros(endo_nbr,h);
K=nan(endo_nbr,pp,h,h);
for snow=1:h
    for slag=1:h
        [a01tt(:,snow,slag),P01tt(:,:,snow,slag),K(:,:,snow,slag)]=kalman_update(a01(:,snow,slag),P01(:,:,snow,slag),...
            iF(:,:,snow,slag),v(:,snow,slag),obs_id);
        att(:,snow)=att(:,snow)+PAI01_tt(snow,slag)*a01tt(:,snow,slag);
    end
    att(:,snow)=att(:,snow)/PAItt(snow);
end
Ptt=zeros(endo_nbr,endo_nbr,h);
for snow=1:h
    for slag=1:h
        Ptt(:,:,snow)=Ptt(:,:,snow)+PAI01_tt(snow,slag)*(P01tt(:,:,snow,slag)+...
            (att(:,snow)-a01tt(:,snow,slag))*transpose(att(:,snow)-a01tt(:,snow,slag)));
    end
    Ptt(:,:,snow)=Ptt(:,:,snow)/PAItt(snow); %symmetrize()
end
end

