function parts=partition_variables(LLI)
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

% this function assumes that LLI is ordered alphabetically
parts=struct();
static = []; pred = []; both = []; frwrd = [];
if ~isempty(LLI)
    static = find(~LLI(:,1) & ~LLI(:,3));
    pred = find(~LLI(:,1) & LLI(:,3));
    both = find(LLI(:,1) & LLI(:,3));
    frwrd = find(LLI(:,1) & ~LLI(:,3));
end
parts.static=static(:)';
parts.pred=pred(:)';
parts.both=both(:)';
parts.frwrd=frwrd(:)';
parts.order_var=[parts.static,parts.pred,parts.both,parts.frwrd];
parts.inv_order_var(parts.order_var)=1:numel(parts.order_var);
end
