function varargout=state_matrices(pairings,varargin)
% examples of calls:
% zeroth_order=state_matrices(pairings,G0)
% [zeroth_order,first_order]=state_matrices(pairings,G0,G1)
% [zeroth_order,first_order,second_order]=state_matrices(pairings,G0,G1,G2)
% 
% pairings is a structure whose number of elements is the order + 1. It has
% fields pairings
% p --> plus or leads: nc %--> f?
% c --> current or 0: nc %--> ok
% l --> lags : nc %--> b?
% x --> shocks : nx %--> e?
% t --> parameters (theta): nt %--> p?

order=length(varargin)-1;
if order<0
    error('at least one jacobian should be entered')
end
if order>2
    error('orders greater than 2 not implemented yet')
end

neqtns=pairings(1).size(1);
nc=pairings(2).pairings.c.ncols;
nx=pairings(2).pairings.x.ncols;
nt=pairings(2).pairings.t.ncols;
h=size(varargin{1},2);

zeroth_order=[];
first_order=[];
second_order=[];
initialize_output();

G0=varargin{1};
if order>0
    G1=varargin{2};
    if order>1
        G2=varargin{3};
    end
end
clear varargin

zeroth_order.G0=G0;
for s0=1:h
    for s1=1:h
        if order>0
            first_order.Gp{s0,s1}(:,real(Index1.p))=G1{s0,s1}(:,imag(Index1.p));
            first_order.Gc{s0,s1}(:,real(Index1.c))=G1{s0,s1}(:,imag(Index1.c));
            first_order.Gl{s0,s1}(:,real(Index1.l))=G1{s0,s1}(:,imag(Index1.l));
            first_order.Gx{s0,s1}(:,real(Index1.x))=G1{s0,s1}(:,imag(Index1.x));
            first_order.Gt{s0,s1}(:,real(Index1.t))=G1{s0,s1}(:,imag(Index1.t));
            if order>1
                second_order.Gpp{s0,s1}(:,real(Index2.pp))=G2{s0,s1}(:,imag(Index2.pp));
                second_order.Gpc{s0,s1}(:,real(Index2.pc))=G2{s0,s1}(:,imag(Index2.pc));
                second_order.Gpl{s0,s1}(:,real(Index2.pl))=G2{s0,s1}(:,imag(Index2.pl));
                second_order.Gpx{s0,s1}(:,real(Index2.px))=G2{s0,s1}(:,imag(Index2.px));
                second_order.Gpt{s0,s1}(:,real(Index2.pt))=G2{s0,s1}(:,imag(Index2.pt));
                second_order.Gcp{s0,s1}(:,real(Index2.cp))=G2{s0,s1}(:,imag(Index2.cp));
                second_order.Gcc{s0,s1}(:,real(Index2.cc))=G2{s0,s1}(:,imag(Index2.cc));
                second_order.Gcl{s0,s1}(:,real(Index2.cl))=G2{s0,s1}(:,imag(Index2.cl));
                second_order.Gcx{s0,s1}(:,real(Index2.cx))=G2{s0,s1}(:,imag(Index2.cx));
                second_order.Gct{s0,s1}(:,real(Index2.ct))=G2{s0,s1}(:,imag(Index2.ct));
                second_order.Glp{s0,s1}(:,real(Index2.lp))=G2{s0,s1}(:,imag(Index2.lp));
                second_order.Glc{s0,s1}(:,real(Index2.lc))=G2{s0,s1}(:,imag(Index2.lc));
                second_order.Gll{s0,s1}(:,real(Index2.ll))=G2{s0,s1}(:,imag(Index2.ll));
                second_order.Glx{s0,s1}(:,real(Index2.lx))=G2{s0,s1}(:,imag(Index2.lx));
                second_order.Glt{s0,s1}(:,real(Index2.lt))=G2{s0,s1}(:,imag(Index2.lt));
                second_order.Gxp{s0,s1}(:,real(Index2.xp))=G2{s0,s1}(:,imag(Index2.xp));
                second_order.Gxc{s0,s1}(:,real(Index2.xc))=G2{s0,s1}(:,imag(Index2.xc));
                second_order.Gxl{s0,s1}(:,real(Index2.xl))=G2{s0,s1}(:,imag(Index2.xl));
                second_order.Gxx{s0,s1}(:,real(Index2.xx))=G2{s0,s1}(:,imag(Index2.xx));
                second_order.Gxt{s0,s1}(:,real(Index2.xt))=G2{s0,s1}(:,imag(Index2.xt));
                second_order.Gtp{s0,s1}(:,real(Index2.tp))=G2{s0,s1}(:,imag(Index2.tp));
                second_order.Gtc{s0,s1}(:,real(Index2.tc))=G2{s0,s1}(:,imag(Index2.tc));
                second_order.Gtl{s0,s1}(:,real(Index2.tl))=G2{s0,s1}(:,imag(Index2.tl));
                second_order.Gtx{s0,s1}(:,real(Index2.tx))=G2{s0,s1}(:,imag(Index2.tx));
                second_order.Gtt{s0,s1}(:,real(Index2.tt))=G2{s0,s1}(:,imag(Index2.tt));
            end
        end
    end
end
varargout=cell(1,order+1);
varargout{1}=zeroth_order; clear zeroth_order
if order>0
    varargout{2}=first_order; clear first_order
    if order>1
        varargout{2}=second_order; clear second_order
    end
end

    function Index=load_pairings(ordr)
        Index=pairings(ordr+1).pairings;
        fields=fieldnames(Index);
        for ifield=1:numel(fields)
            ff=fields{ifield};
            Index.(ff)=Index.(ff).pairs;
        end
    end

    function initialize_output()
        zeroth_order=struct('G0',{repmat({sparse(neqtns,1)},h,h)});
        first_order=[];
        second_order=[];
        if order>0
            Index1=load_pairings(1);
            first_order=struct(...
                'Gp',{repmat({sparse(neqtns,nc)},h,h)},...
                'Gc',{repmat({sparse(neqtns,nc)},h,h)},...
                'Gl',{repmat({sparse(neqtns,nc)},h,h)},...
                'Gx',{repmat({sparse(neqtns,nx)},h,h)},...
                'Gt',{repmat({sparse(neqtns,nt)},h,h)}...
                );
            if order>1
                Index2=load_pairings(2);
                second_order=struct(...
                    'Gpp',{repmat({sparse(neqtns,nc*nc)},h,h)},...
                    'Gpc',{repmat({sparse(neqtns,nc*nc)},h,h)},...
                    'Gpl',{repmat({sparse(neqtns,nc*nc)},h,h)},...
                    'Gpx',{repmat({sparse(neqtns,nc*nx)},h,h)},...
                    'Gpt',{repmat({sparse(neqtns,nc*nt)},h,h)},...
                    'Gcp',{repmat({sparse(neqtns,nc*nc)},h,h)},...
                    'Gcc',{repmat({sparse(neqtns,nc*nc)},h,h)},...
                    'Gcl',{repmat({sparse(neqtns,nc*nc)},h,h)},...
                    'Gcx',{repmat({sparse(neqtns,nc*nx)},h,h)},...
                    'Gct',{repmat({sparse(neqtns,nc*nt)},h,h)},...
                    'Glp',{repmat({sparse(neqtns,nc*nc)},h,h)},...
                    'Glc',{repmat({sparse(neqtns,nc*nc)},h,h)},...
                    'Gll',{repmat({sparse(neqtns,nc*nc)},h,h)},...
                    'Glx',{repmat({sparse(neqtns,nc*nx)},h,h)},...
                    'Glt',{repmat({sparse(neqtns,nc*nt)},h,h)},...
                    'Gxp',{repmat({sparse(neqtns,nx*nc)},h,h)},...
                    'Gxc',{repmat({sparse(neqtns,nx*nc)},h,h)},...
                    'Gxl',{repmat({sparse(neqtns,nx*nc)},h,h)},...
                    'Gxx',{repmat({sparse(neqtns,nx*nx)},h,h)},...
                    'Gxt',{repmat({sparse(neqtns,nx*nt)},h,h)},...
                    'Gtp',{repmat({sparse(neqtns,nt*nc)},h,h)},...
                    'Gtc',{repmat({sparse(neqtns,nt*nc)},h,h)},...
                    'Gtl',{repmat({sparse(neqtns,nt*nc)},h,h)},...
                    'Gtx',{repmat({sparse(neqtns,nt*nx)},h,h)},...
                    'Gtt',{repmat({sparse(neqtns,nt*nt)},h,h)} ...
                    );
            end
        end
    end
end
