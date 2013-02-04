function [T,R,SS]=recast_loose_commitment_solution_into_markov_switching(obj,T,R,SS)
% this function puts the loose commitment solution into a markov switching
% form for estimation and simulation.

error(nargchk(4,4,nargin,'struct'))

% should write a function called  that does
% what I do here, so that I can be used for simulation and for the
% calculation of moments too.
if obj.is_optimal_policy_model && obj.NumberOfRegimes==2
    if obj.options.lc_reconvexify
        error([mfilename,':: reconvexification is not compatible with multiple regimes'])
    end
    T=repmat(T,[1,1,2]);
    % here is the switch: in the alternative regime, we ignore past
    % promises
    T(:,obj.is_lagrange_multiplier,2)=0;
    R=repmat(R,[1,1,1,2]);
    SS=SS(:,ones(1,2));
end

% % function [T,R,SS,Q,H]=recast_loose_commitment_solution_into_markov_switching(obj,T,R,SS,Q,H)
% % % this function puts the loose commitment solution into a markov switching
% % % form for estimation and simulation.
% %
% % error(nargchk(4,6,nargin,'struct'))
% %
% % if nargin<6
% %     H=[];
% %     if nargin<5
% %         Q=1;
% %     end
% % end
% %
% % % should write a function called  that does
% % % what I do here, so that I can be used for simulation and for the
% % % calculation of moments too.
% % if obj.is_optimal_policy_model && ~obj.options.lc_reconvexify
% %     gam=obj.planner_commitment;
% %     if ~ismember(gam,[0,1])
% %         T=repmat(T,[1,1,2]);
% %         % here is the switch: in the alternative regime, we ignore past
% %         % promises
% %         T(:,obj.is_lagrange_multiplier,2)=0;
% %         R=repmat(R,[1,1,1,2]);
% %         SS=SS(:,ones(1,2));
% %         if ~isempty(H)
% %             H=repmat(H,[1,1,2]);
% %         end
% %         if obj.NumberOfRegimes==1
% %             Q=[gam,1-gam
% %                 1,0];
% %         else
% %             % well, Q must be taken care of. All we need to do is to
% %             % check that the number of regimes does not exceed 2
% %             assert(obj.NumberOfRegimes==2,'In loose commitment, the number of regimes cannot exceed 2')
% %         end
% %     end
% % end
