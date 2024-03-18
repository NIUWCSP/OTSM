function [Gn_block_matrix,Tn_block_matrix,zn_block_vector,H_t_f]=Generate_Matched_Filter_GS_matrices(N,M,G,r,ChanEst)

ChanEstGrid = repmat(ChanEst,N/size(ChanEst,1),M/size(ChanEst,2));
% Generate time-domain GS iteration matrices for low complexity iterative detection
Gn_block_matrix=zeros(M,M,N);
Tn_block_matrix=zeros(M,M,N);
Qn_block_matrix=zeros(M,M,N);
zn_block_vector=zeros(M,N);

H_t_f=zeros(N,M); % Time-frequency single tap channel matrix
Fn=dftmtx(M);%离散傅里叶变换矩阵
Fn=Fn./norm(Fn);
for n=1:N
    rn=r((n-1)*M+1:n*M);
    Gn_block_matrix(:,:,n)=G((n-1)*M+1:n*M,(n-1)*M+1:n*M);
    %Gn=Gn_block_matrix(:,:,n)./ChanEstGrid;
    Gn=Gn_block_matrix(:,:,n);
    H_t_f(n,1:M)=diag(Fn*Gn*Fn').';  % Generate time-frequency channel matrix for low complexity initial estimate using equation (20) in [R2]
    Rn=Gn'*Gn;    %2.22
    Dn=diag(diag(Rn));
    Ln=tril(Rn,-1);
    Un=triu(Rn,1);
    Qn=(Dn+Ln)^(-1);
    Tn=Qn*Un;
    Tn_block_matrix(:,:,n)=Tn;
    Qn_block_matrix(:,:,n)=Qn;        
    zn_block_vector(:,n)=Qn*Gn'*rn;
end
end