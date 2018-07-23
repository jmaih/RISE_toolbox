function [myirfs]=irf(self,shock_names,irf_periods,params,Rfunc)
% Compute impulse response function from the given parameter values
%
% ::
%
%    myirfs = irf(self, shock_names, irf_periods, params, Rfunc);
%
% Args:
%    self (var object): var object
%    shock_names (cellstr): shocks to compute IRFs (has to be consistent with the names given in :func:`identification <var.identification>`)
%    irf_periods (int): number of periods to compute IRFs (default: 40)
%    params (vector): parameter values of the model. If empty MLE/posterior-mode values are used.
%    Rfunc (function handle): identification function. This is an output of :func:`identification <var.identification>`. (default: choleski identification)
%
% Returns:
%    : struct containing IRFs
%
% Note:
%    Only successful IRFs are returned. If the structure does not return some IRFs make sure that IRFs properly identified.
%

n=nargin;

set_defaults()

[myirfs,info]=vartools.irf(self,irf_periods,params,Rfunc); %#ok<ASGLU>

tmp=struct();

shock_names=abstvar.create_variable_names(self.nvars,'shock',shock_names);

for ishock=1:self.nvars

    sname=shock_names{ishock};

    batchShock=(ishock-1)*self.ng+1:ishock*self.ng;

    for iv=1:self.nvars

        vname=self.endogenous{iv};

        batchVars=(iv-1)*self.ng+1:iv*self.ng;

        for ishock2=1:numel(batchShock)

            possh=batchShock(ishock2);

            if self.is_panel

            C1=self.members{ishock2};

            end

            for iv2=1:numel(batchVars)

                posv=batchVars(iv2);

                if self.is_panel

                C2=self.members{iv2};

                end

                data=permute(myirfs(posv,:,possh,:,:),[2,5,4,1,3]);

                data=squeeze(data);

                if ishock==1 && iv==1 && ishock2==1 && iv2==1

                    if self.nregs==1

                        regimeNames=[];

                    else

                        regimeNames=abstvar.create_variable_names(self.nregs,'regime');

                    end

                    proto=ts(1,data,regimeNames);

                end

                if self.is_panel

                    tmp.(sname).(C1).(vname).(C2)=reset_data(proto,data,...
                        regimeNames);

                else

                    tmp.(sname).(vname)=reset_data(proto,data,regimeNames);

                end

            end

        end

    end

end

myirfs=tmp;

    function set_defaults()

        if n < 5

            Rfunc=[];

            if n < 4

                params=[];

                if n< 3

                    irf_periods=[];

                    if n<2

                        shock_names=[];

                    end

                end

            end

        end

        if isempty(Rfunc),Rfunc=identification(self,'choleski'); end

        params=solve(self,params);

        if isempty(irf_periods),irf_periods=40; end

    end

end
