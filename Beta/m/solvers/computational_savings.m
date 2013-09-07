function [AA0,AAminus,AAplus,BB,statcols,q]=computational_savings(A0,Aminus,Aplus,enforce,varargin)
% this function separates static variables from dynamic ones. It places all
% the static variables at the beginning and the dynamic ones following
% after.
% statcols is a logical vector indicating the original position of the
% static variables.
% the matrices are transformed in such a way that the dynamic equations
% appear at the bottom and can be solved independently. Once the solution
% for the dynamic variables is found, the solution for the static variables
% can be computed as a function of the dynamic variables.
if nargin<4
    enforce=[];
end
statcols=~any(Aplus) & ~any(Aminus);
% othercols=1-statcols;

if size(statcols,3)>1
    statcols=any(transpose(squeeze(statcols)));
end
% the variables enforced are not statcols
statcols(enforce)=false;

BB=varargin;
if any(statcols)
    [AA0,AAminus,AAplus]=re_order(statcols,A0,Aminus,Aplus);
    [nrows,~,nreg]=size(AA0);
    q=zeros(nrows,nrows,nreg);
    for ii=1:nreg
        [q(:,:,ii),r]=qr(AA0(:,:,ii)); %#ok<NASGU>
        AA0(:,:,ii)=q(:,:,ii)'*AA0(:,:,ii);
        AAminus(:,:,ii)=q(:,:,ii)'*AAminus(:,:,ii);
        AAplus(:,:,ii)=q(:,:,ii)'*AAplus(:,:,ii);
        for jj=1:numel(BB)
            BB{jj}(:,:,ii)=q(:,:,ii)'*BB{jj}(:,:,ii);
        end
    end
else
    [AA0,AAminus,AAplus]=deal(A0,Aminus,Aplus);
    q=1;
end
function varargout=re_order(statcols,varargin)

for ii=1:numel(varargin)
    varargout{ii}=[varargin{ii}(:,statcols,:),varargin{ii}(:,~statcols,:)];
end