function [est_bits,ite,x_data] = Matched_Filter_GS_detector(N,M,M_mod,no,data_grid,Y,H_tf,n_ite,omega,Tn_block_matrix,Gn_block_matrix,zn_block_vector,r,Wn,decision)
%% Initial assignments
%Number of symbols per frame
N_syms_perfram=sum(sum((data_grid>0)));
%Arranging the delay-Doppler grid symbols into an array
data_array=reshape(data_grid,1,N*M);
%finding position of data symbols in the array
[~,data_index]=find(data_array>0);
M_bits=log2(M_mod);
N_bits_perfram = N_syms_perfram*M_bits;
%% initial time-frequency low complexity estimate assuming ideal pulses
if(1)
    Y_tf=fft((Y*Wn)).'; % ISFFT                                              %% equation (18) in [R2]
    X_tf=conj(H_tf).*Y_tf./(H_tf.*conj(H_tf)+no); % single tap equalizer     %% equation (19) in [R2]
    X_est = ifft(X_tf.')*Wn; % SFFT                                          %% equation (21) in [R2]
    X_est=qammod(qamdemod(X_est,M_mod,'gray'),M_mod,'gray');
    X_est=X_est.*data_grid;
    X_tilda_est=X_est*Wn;
end
X_tilda_est=X_tilda_est.*data_grid;

%% Matched Filter Gauss Siedel algorithm
error=zeros(n_ite);
Y_tilda_est=zeros(M,N);
Y_tilda=reshape(r,M,N);
x_soft=zeros(M,N);

for ite=1:n_ite
    for i=1:N
        Tn=Tn_block_matrix(:,:,i);
        zn=zn_block_vector(:,i);
        sn_prev=X_tilda_est(:,i);
        sn_next=-Tn*sn_prev+zn;     %equation (47) in [R1] (or equivalently equation (27) in [R2])
        x_soft(:,i)=sn_next.*data_grid(:,i);
        X_tilda_est(:,i)=(x_soft(:,i));        
    end
    if(decision==1)
        x_m=x_soft*Wn;
        X_tilda_est=(1-omega)*X_tilda_est+omega*((qammod(qamdemod(x_m,M_mod,'gray'),M_mod,'gray').*data_grid)*Wn);  %equation (50) in [R1] (or equivalently equation (27) in [R2])
    end    
    for i=1:N
        Gn=Gn_block_matrix(:,:,i);
        Y_tilda_est(1:M,i)=Gn*(X_tilda_est(1:M,i));
    end
    error(ite)=sum(sum(abs(Y_tilda_est-Y_tilda).^2))./N_syms_perfram;
    if(ite>1)
        if(error(ite)>=error(ite-1))
            break;
        end
    end    
end
if(n_ite==0)
    ite=0;
end

%% detector output likelihood calculations for turbo decode
X_est=X_tilda_est*Wn;
x_est=reshape(X_est,1,N*M);
x_data=x_est(data_index);
est_bits=reshape(qamdemod(x_data,M_mod,'gray','OutputType','bit'),N_bits_perfram,1);

%畫圖
subplot(221);
plot(x_m,'.');title('MFGS detector');axis equal;
end