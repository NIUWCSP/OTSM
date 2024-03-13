function  [G,r,Y] = Receiver(N,M,car_fre,delta_f,T,iesn0,sigma_2,z,Wn)
%NumFFT = 64;%V3 FFT轉換的點數
%NumSyncPreamble = 32;%V3 同步的前綴，Preamble：防干擾+同步+通道估測(已知的頻域資料)
%NumCP = 16;%V3 CP：循環前綴，CP：避免ISI(多路徑干擾)(未知的時域訊號)
max_speed=500;  % km/hr
[chan_coef,delay_taps,Doppler_taps,taps]=Generate_delay_Doppler_channel_parameters(N,M,car_fre,delta_f,T,max_speed);
           
        %% channel output%%%%%
        [G,gs]=Gen_time_domain_channel(N,M,taps,delay_taps,Doppler_taps,chan_coef);
        
        r=zeros(N*M,1);
        noise = sqrt(sigma_2(iesn0)/2) * (randn(size(z)) + 1i*randn(size(z)));
        l_max=max(delay_taps);
        for q=0:N*M-1
            for l=0:l_max
                if(q>=l)
                    r(q+1)=r(q+1)+gs(l+1,q+1)*z(q-l+1);  %equation (24) in [R1]
                end
            end
        end
        r=r+noise;
           
        %% OTSM demodulation%%%%
        Y_tilda=reshape(r,M,N);     %equation (11) in [R1]
        Y = Y_tilda*Wn;             %equation (12) in [R1]
        
        










end