function [chan_coef,delay_taps,Doppler_taps,taps]=Generate_delay_Doppler_channel_parameters(N,M,car_fre,delta_f,T,max_speed)
one_delay_tap = 1/(M*delta_f);
one_doppler_tap = 1/(N*T);

%delays = [0 30 70 90 110 190 410]*10^(-9);%EPA model 
delays = [0 30 150 310 370 710 1090 1730 2510]*10^(-9);%EVA model
% delays = [0 50 120 200 230 500 1600 2300 5000]*10^(-9);%ETU model


taps = length(delays);% number of delay taps
delay_taps = round(delays/one_delay_tap);%assuming no fraction for the delay

%pdp = [0 -1.0 -2.0 -3.0 -8.0 -17.2 -20.8];% EPA power delay profile
pdp = [0 -1.5 -1.4 -3.6 -0.6 -9.1 -7.0 -12.0 -16.9];% EVA power delay profile
%pdp = [-1 -1 -1 0 0 0  -3 -5 -7];% ETU power delay profile


pow_prof = 10.^(pdp/10);
pow_prof = pow_prof/sum(pow_prof);%normalization of power delay profile
chan_coef = sqrt(pow_prof).*(sqrt(1/2) * (randn(1,taps)+1i*randn(1,taps)));%channel coef. for each path
max_UE_speed = max_speed*(1000/3600);
Doppler_vel = (max_UE_speed*car_fre)/(299792458);
max_Doppler_tap = Doppler_vel/one_doppler_tap;
Doppler_taps = (max_Doppler_tap*cos(2*pi*rand(1,taps)));%Doppler taps using Jake's spectrum
end