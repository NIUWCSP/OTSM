function [ip, upsample] = pluto()

clearvars -except times;close all;warning off; %預設環境
set(0,'defaultfigurecolor','w'); 
%加入path
addpath ..\..\library 
addpath ..\..\library\matlab 
addpath ..\..\code\matlab\OFDM

%刪除.mat
if(0)
    Delete_mat;
end

% 載入資料
if exist('ScattorData.mat','file')
    load('ScattorData.mat');
else
    AllRxDataSymbEqAverage = zeros(0,0); %儲存Scattor資料
end
if exist('BERData.mat','file')
    load('BERData.mat');
else
    AllBERData = zeros(0,0); %儲存BER資料
end
%設定pluto IP
ip = '192.168.2.1';

%設定與進入TX函式
upsample=4; %過取樣取4倍，數位還原類比後比較可以不失真
end