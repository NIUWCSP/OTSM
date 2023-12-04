function [G,gs]=Gen_time_domain_channel(N,M,P,delay_taps,Doppler_taps,chan_coef)
z=exp(1i*2*pi/N/M);
l_max=max(delay_taps);
gs=zeros(l_max+1,N*M);      
G=zeros(N*M,N*M);

for q=0:N*M-1
    for i=1:P
        g_i=chan_coef(i);
        l_i=delay_taps(i);
        k_i=Doppler_taps(i);        
        gs(l_i+1,q+1)=gs(l_i+1,q+1)+g_i*z^(k_i*(q-l_i));  % equation (22) in [R1]
    end    
end

for q=0:N*M-1
    for l=0:l_max
        if(q>=l)
            G(q+1,q+1-l)=gs(l+1,q+1);
        end
    end
end
end