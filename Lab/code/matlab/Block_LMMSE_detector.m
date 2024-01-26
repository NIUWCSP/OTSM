function [est_bits,x_data] = Block_LMMSE_detector(N,M,M_mod,noise_var,data_grid,Gn_block_matrix,r,Wn)
%% Initial assignments
%Number of symbols per frame
N_syms_perfram=sum(sum((data_grid>0)));
%Arranging the delay-Doppler grid symbols into an array
data_array=reshape(data_grid,1,N*M);
%finding position of data symbols in the array
[~,data_index]=find(data_array>0);
M_bits=log2(M_mod);
N_bits_perfram = N_syms_perfram*M_bits;
sn_block_est=zeros(M,N);
%% Gauss Siedel SOR
for i=1:N    
    rn=r((i-1)*M+1:i*M);
    Gn=Gn_block_matrix(:,:,i);
    Rn=Gn'*Gn;%共變異數矩陣
    sn_block_est(:,i)=(Rn+noise_var.*eye(M))^(-1)*(Gn'*rn);
end
X_tilda_est=sn_block_est;
%% detector output
X_est=X_tilda_est*Wn;
x_est=reshape(X_est,1,N*M);
x_data=x_est(data_index);
est_bits=reshape(qamdemod(x_data,M_mod,'gray','OutputType','bit'),N_bits_perfram,1);
end