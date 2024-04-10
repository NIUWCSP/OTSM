function [est_bits,x_data] = TF_single_tap_equalizer(N,M,M_mod,noise_var,data_grid,Y,H_tf,Wn)
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
    Y_tf=fft((Y*Wn)).'; % ISFFT
    X_tf=conj(H_tf).*Y_tf./(H_tf.*conj(H_tf)+noise_var); % single tap equalizer
    X_est = ifft(X_tf.')*Wn; % SFFT
    X_est=qammod(qamdemod(X_est,M_mod,'gray'),M_mod,'gray');
    X_est=X_est.*data_grid;
    X_tilda_est=X_est*Wn;
end
X_tilda_est=X_tilda_est.*data_grid;

%% detector output likelihood calculations for turbo decode
X_est=X_tilda_est*Wn;
x_est=reshape(X_est,1,N*M);
x_data=x_est(data_index);
est_bits=reshape(qamdemod(x_data,M_mod,'gray','OutputType','bit'),N_bits_perfram,1);

%畫圖
subplot(222);
plot(ifft(X_tf.')*Wn,'.');title('single tap equalizer');axis equal;
end