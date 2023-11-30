%% MK Mutation Test_Main
%Yanpeng Yang, 2023.10.10
%yangyanpeng@nieer.ac.cn

function [UF,UB] = MKbreak(time_series)
n = length(time_series);
%%  ---------------------------------Positive sequence computation--------------------------------
% Defining Variable UF£¬Lenth = n£¬Initial value = 0£»
UF=zeros(size(time_series));

%E = n*(n-1)/4;                         % E is the mean of the cumulative number s(k)
%Var = n*(n-1)*(2*n+5)/72;              % Var is the variance of the cumulative number s(k)

% To define the order column, r(i) records the number of values at the i th moment 
% whose value is greater than the number at the j moment (where j<=i)
for i= 1:n
    r1(i) = sum(time_series(i)>time_series(1:i));
end

% Define the sequence of cumulators s, s(k) records the i th moment (where i<=k), 
% whose value is greater than the sum of the number of values at the moment j (where j<=i)
s = zeros(size(time_series));
% k starts at 2, because according to the statistic UF(k) formula, when k=1, s(1), E(1) and Var(1) are all 0, 
% and UF is meaningless at this time, so UFk(1) is set to 0 in the formula
for k = 2:n
   s(k) = sum(r1(1:k));

   E = k*(k-1)/4;                     % The mean of s(k)
   Var = k*(k-1)*(2*k+5)/72;          % The variance of s(k)
   UF(k) = (s(k)-E)/sqrt(Var);
end

%% ---------------------------------Inverse column computation--------------------------------
% Reverse the sample by time series, constructing the reverse sequence time_series2
for i=1:n
    time_series2(i)=time_series(n-i+1);
end
% Define the reverse statistic UB with length =n and initial value =0
UB = zeros(size(time_series2));
for i= 1:n
    r2(i) = sum(time_series2(i)>time_series2(1:i));
end

% Define the reverse cumulative sequence s2 with length =n and initial value =0
s2 = zeros(size(time_series2));
% i starts from 2, because according to the statistic UB(k) formula, when i=1, s2(1), E(1) and Var(1) are all 0, 
% and UB is meaningless at this time, so UB(1)=0 is set in the formula
for k = 2:n
   s2(k) = sum(r2(1:k));
   
   E = k*(k-1)/4;                     % The mean of s2(k)
   Var = k*(k-1)*(2*k+5)/72;          % Variance of s2(k)
% Since the accumulation method is still used in the construction of the accumulative quantity S2 of the reverse sequence, that is, r adds 1 when the latter is greater than the former,
% The size of r represents the size of an upward trend, and after the sequence is reversed, it should show the trend opposite to the original sequence.
% Therefore, the statistical formula (S(i)-E(i))/sqrt(Var(i)) should not be changed when using the summative method to count S2 series.
% But the statistic UBk should be negative to represent the trend of the correct reverse sequence
  UB(k) = -(s2(k)-E)/sqrt(Var);
end

%% ---------------------------------Plot--------------------------------------------------
% At this time, the trend statistics of the reverse order list on the reverse order time shown in 
% the previous step to UB is the same as that of UF when plotting to find abrupt points, the two 
% curves should have the same time axis, so the result statistics UB is reversed according to the 
% time series, and the UB of the positive time order is obtained for plotting

for i=1:n
   UB2(i)=UB(n-i+1);
end
x = 1:n;

plot(x,UF,'-','color',[0.18,0.55,0.34],'linewidth',1.5);
hold on
plot(x,UB2,'m--','linewidth',1.5);
%plot(x,UB,'m-','linewidth',1.5);

plot(x,2.58*ones(n,1),'-','color',[0.74,0.71,0.42],'linewidth',1);
plot(x,1.96*ones(n,1),'-.','color',[0.74,0.71,0.42],'linewidth',1);
%plot(x,1.28*ones(n,1),':','color',[0.74,0.71,0.42],'linewidth',1);
plot(x,0*ones(n,1),'k-','linewidth',0.6);
%plot(x,-1.28*ones(n,1),':','color',[0.74,0.71,0.42],'linewidth',1);
plot(x,-1.96*ones(n,1),'-.','color',[0.74,0.71,0.42],'linewidth',1);
plot(x,-2.58*ones(n,1),'-','color',[0.74,0.71,0.42],'linewidth',1);
grid(gca,'minor')

axis([min(x),max(x),-8,8]);
legend('UF','UB','¦Á=0.01','¦Á=0.05');
set(gca,'Fontsize',12)
set(gca,'ytick',-8:2:8)
xlabel('{\itt} (year)','FontName','TimesNewRoman','FontSize',15);
ylabel('Statistic variables','FontName','TimesNewRoman','Fontsize',15);

end

