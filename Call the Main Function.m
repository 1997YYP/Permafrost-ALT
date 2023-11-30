%% MK Mutation Test_Call main function
%Yanpeng Yang, 2023.10.10
%yangyanpeng@nieer.ac.cn

clc;
time_series = xlsread('H:\040 Permafrost Degradation\016 Consistency Analysis\1981-2021_Consistency Analysis.xlsx','NDVI Anomaly','AP288:BO288');
A=MKbreak(time_series);
disp(A);