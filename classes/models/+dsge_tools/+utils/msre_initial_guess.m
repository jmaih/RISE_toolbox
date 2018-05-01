function T0=msre_initial_guess(A0,dpb_minus,dbf_plus,solve_initialization,old_T0)
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

if nargin==0
    T0=struct();
    return
end
n=size(A0{1},1);
h=numel(A0);
npb=size(dpb_minus{1},2);
T0=zeros(n,npb,h);
if ~isempty(old_T0)
    for ireg=1:h
        T0(:,:,ireg)=old_T0{ireg}(:,1:npb);
    end
    return
end

switch solve_initialization
    case {'zeros'}
    case {'backward'}
        %     warnstate=warning('query','all');
        %     warning('off','MATLAB:singularMatrix')
        for i_state=1:h
            % the solution that corresponds to the backward-looking model [default]
            % this should accelerate the solving and give more accurate results since
            % we strictly don't need to solve for the non-state columns i.e. the columns
            % where Aminus is zero
%             T0(:,:,i_state)=-A0{i_state}\dpb_minus{i_state};
            T0(:,:,i_state)=-pinv(full(A0{i_state}))*dpb_minus{i_state};
%             if any(any(isnan(T0(:,:,i_state))))
%                 T0(:,:,i_state)=-pinv(full(A0{i_state}))*dpb_minus{i_state};
%                 % the A0 matrix can be singular, for instance in the
%                 % case of the zlb
%             end
        end
        %     warning(warnstate)
    case {'random'}
        level=3;
        switch level
            case 0
                T0=randn(n,npb,h);
            case 1
                sn=solvent_norm();
                AT=sn^2*randn(n,npb,h);
                for istate=1:h
                    QAT=AT(:,:,jstate)+A0{istate};
                    T0(:,:,istate)=-QAT\dpb_minus{istate};
                end
            case 2
                Tbkw=msre_initial_guess('back_init');
                T0=Tbkw.*rand(n,npb,h);
            case 3
                sn=solvent_norm();
                T0=sn*randn(n,npb,h);
        end
end
    function sn=solvent_norm()
        n_a0=0;n_aplus=0;n_aminus=0;
        for ii_state=1:h
            n_a0=max(n_a0,norm(full(A0{ii_state})));
            aplus=0;
            for s1=1:h
                aplus=aplus+dbf_plus{ii_state,s1};
            end
            n_aplus=max(n_aplus,norm(aplus));
            n_aminus=max(n_aminus,norm(full(dpb_minus{ii_state})));
        end
        sn=(n_a0+sqrt(n_a0^2+4*n_aplus*n_aminus))/(2*n_aplus);
    end
end
