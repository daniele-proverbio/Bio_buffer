%% Plotting final results


clear 
close all
clc

%% Loading
% Load the statistics indicators calculated in analysis.m
% statistics contains 

n2 = load('statistics2.mat');
n3 = load('statistics3.mat');
n4 = load('statistics4.mat');
n5 = load('statistics5.mat');
n8 = load('statistics8.mat');

%%
cc0 = flip(0.15:0.005:0.35);
stop=40;

%% Variance

figure()
hold on
plot(cc0(1:stop),n2.statistics.mean_var(1:stop),linewidth=1.2)
plot(cc0(1:stop),n3.statistics.mean_var(1:stop),linewidth=1.2)
plot(cc0(1:stop),n4.statistics.mean_var(1:stop),linewidth=1.2)
plot(cc0(1:stop),n5.statistics.mean_var(1:stop),linewidth=1.2)
plot(cc0(1:stop),n8.statistics.mean_var(1:stop),linewidth=1.2)
ylabel('Var',fontsize=20,Interpreter='latex',fontweight='bold')
xlabel('$c-c_0$',fontsize=20,Interpreter='latex',fontweight='bold')
ax = gca;
ax.FontSize = 18;
%legend({'n=2','n=3','n=4','n=5','n=8'},FontSize=14,Location='northeast')

%% Autocorrelation

figure()
hold on
plot(cc0(1:stop),n2.statistics.mean_AC(1:stop),linewidth=1.2)
plot(cc0(1:stop),n3.statistics.mean_AC(1:stop),linewidth=1.2)
plot(cc0(1:stop),n4.statistics.mean_AC(1:stop),linewidth=1.2)
plot(cc0(1:stop),n5.statistics.mean_AC(1:stop),linewidth=1.2)
plot(cc0(1:stop),n8.statistics.mean_AC(1:stop),linewidth=1.2)
ylabel('AC(1)',fontsize=20,Interpreter='latex',fontweight='bold')
xlabel('$c-c_0$',fontsize=20,Interpreter='latex',fontweight='bold')
ax = gca;
ax.FontSize = 18;
%legend({'n=2','n=3','n=4','n=5','n=8'},FontSize=14,Location='southwest')



%% p-value var
N= 5;
p_value_thresholds = [0.01, 0.05, 0.1];

mean_p_value = (n2.statistics.p_value + n3.statistics.p_value + n4.statistics.p_value + n5.statistics.p_value + n8.statistics.p_value)/N;
matrix_p_value = [n2.statistics.p_value; n3.statistics.p_value; n4.statistics.p_value; n5.statistics.p_value; n8.statistics.p_value];
stderr_pvalue = std(matrix_p_value)/sqrt(N);

figure()

x1 = subplot(2,1,1);
hold on
plot(cc0(2:stop), mean_p_value(2:end),'k',linewidth=1.5)
plot(cc0(2:stop), mean_p_value(2:end) + stderr_pvalue(2:end),'-.k',linewidth=1.2)
plot(cc0(2:stop), mean_p_value(2:end) - stderr_pvalue(2:end),'-.k',linewidth=1.2)
ax = gca;
ax.FontSize = 18;
set(gca,'XTickLabel',[]);
%ylabel('p-value on Var',fontsize=30,Interpreter='latex')
%xlabel('$c-c_0$',fontsize=30,Interpreter='latex')
yline(0,'-',linewidth=1.2)
yline(p_value_thresholds(2),'--',linewidth=1.2)
yline(p_value_thresholds(3),'--',linewidth=1.2)
legend('Average, over Var','Std err','','','',FontSize=20,Location='northwest')

p_v_var.cc0 = cc0(2:stop);
p_v_var.mean = mean_p_value(2:end);
p_v_var.stderr = stderr_pvalue(2:end);
save(['p_v_var.mat'],'p_v_var', '-v7.3')

% p-value AC
N= 5;
mean_p_valueAC = (n2.statistics.p_valueAC(1:stop) + n3.statistics.p_valueAC(1:stop) + n4.statistics.p_valueAC(1:stop) + n5.statistics.p_valueAC(1:stop) + n8.statistics.p_valueAC(1:stop))/N;
matrix_p_valueAC = [n2.statistics.p_valueAC(1:stop); n3.statistics.p_valueAC(1:stop); n4.statistics.p_valueAC(1:stop); n5.statistics.p_valueAC(1:stop); n8.statistics.p_valueAC(1:stop)];
stderr_pvalueAC = std(matrix_p_valueAC)/sqrt(N);

x2 = subplot(2,1,2);
hold on
plot(cc0(2:stop), mean_p_valueAC(2:end),'k',linewidth=1.5)
plot(cc0(2:stop), mean_p_valueAC(2:end) + stderr_pvalueAC(2:end),'-.k',linewidth=1.2)
plot(cc0(2:stop), mean_p_valueAC(2:end) - stderr_pvalueAC(2:end),'-.k',linewidth=1.2)
ax = gca;
ax.FontSize = 18;
ylabel('p-value',fontsize=30,Interpreter='latex')
xlabel('$c-c_0$',fontsize=30,Interpreter='latex')
yline(p_value_thresholds(2),'--',linewidth=1.2)
yline(p_value_thresholds(3),'--',linewidth=1.2)
legend('Average, over AC(1)','Std err','','','','',FontSize=20,Location='northwest')

p_v_ac.cc0 = cc0(2:stop);
p_v_ac.mean = mean_p_valueAC(2:end);
p_v_ac.stderr = stderr_pvalueAC(2:end);
save(['p_v_ac.mat'],'p_v_ac', '-v7.3')

p1 = get(x1, 'Position');
p2 = get(x2, 'Position');
p1(2) = p2(2)+p2(4);
set(x1, 'pos', p1);
xlabel('$c-c_0$',fontsize=30,Interpreter='latex')