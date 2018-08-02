function [y,newp,retcode]=steady_state_4_Canonical_Const(obj,y,p,d,id) %#ok<INUSL>
% steady_state_4_Canonical_Const --  computes the steady state of ... analytically
%
% ::
%
%   [y,newp,retcode]=steady_state_4_Canonical_Const(obj,y,p,d,id)
%
% Args:
%    obj (rise | dsge): model object (not always needed)
%    y (vector): endo_nbr x 1 vector of initial steady state
%    pp (struct): parameter structure
%    d (struct): definitions
%    id (vector): location of the variables to calculate
%
% Returns:
%    :
%
%    - **y** []: endo_nbr x 1 vector of updated steady state
%    - **newp** [struct]: structure containing updated parameters if any
%    - **retcode** [0|number]: return 0 if there are no problems, else return
%      any number different from 0
%
% Note:
%
%    - this is new approach has three main advantages relative to the previous
%      one:
%
%      - The file is valid whether we have many regimes or not
%      - The user does not need to know what regime is being computed
%      - It is in sync with the steady state model
%

retcode=0;
if nargin==1
    % list of endogenous variables to be calculated
    %----------------------------------------------
    y={'DPQ_P_NW','Y','ZGDP','ZI','ZPAI','ZY','D_GDP_NW','I','PAI','R','RN3M_NW'};
    % list of parameters to be computed during steady state calculation
    %-------------------------------------------------------------------
    newp={};
else
    % no parameters to update or create in the steady state file
    %-----------------------------------------------------------
    newp=[];

    % computation of the steady state
    %--------------------------------
    DPQ_P_NW=p.paiss;
    Y=0;
    ZGDP=p.gyss;
    ZI=0;
    ZPAI=0;
    ZY=0;
    D_GDP_NW=p.gyss;
    I=0;
    PAI=0;
    R=0;
    RN3M_NW=p.iss;

    ys=[DPQ_P_NW , Y , ZGDP , ZI , ZPAI , ZY , D_GDP_NW , I , PAI , R , RN3M_NW ]' ;

    % check the validity of the calculations
    %----------------------------------------
    if ~utils.error.valid(ys)
        retcode=1;
    else
        % push the calculations
        %----------------------
        y(id)=ys;
    end
end

end
