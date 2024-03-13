function [z] = Transmitter(N,M,M_mod,M_bits,data_grid,N_syms_perfram,Wn)
%NumFFT = 64;%V3 FFT轉換的點數
%NumSyncPreamble = 32;%V3 同步的前綴，Preamble：防干擾+同步+通道估測(已知的頻域資料)
%NumCP = 16;%V3 CP：循環前綴，CP：避免ISI(多路徑干擾)(未知的時域訊號)
%% random input bits generation%%%%%
        trans_info_bit = randi([0,1],N_syms_perfram*M_bits,1);%trans_info_bit =TxDataBits
        %%2D QAM symbols generation %%%%%%%%
        data=qammod(reshape(trans_info_bit,M_bits,N_syms_perfram), M_mod,'gray','InputType','bit');%data=1*3840        
        %data=qammod(reshape(trans_info_bit,M_bits,N_syms_perfram), M_mod, 0,'gray','bit');  data=2*3840      
        X = Generate_2D_data_grid(N,M,data,data_grid);
        
        
        %% OTSM modulation%%%%
        X_tilda=X*Wn;               %equation (6) in [R1]
        z = reshape(X_tilda,N*M,1); %equation (7) in [R1]