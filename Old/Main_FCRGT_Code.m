%% Main Code for Functional Cointegrating Regression between GRP growth rate and Temperature Anomaly(FCRGT)
%
%                                                            September 2025
%                                            Kyungsik Nam(ksnam@hufs.ac.kr)

clear; clc; close all; 

cd('C:\Users\noran\Dropbox\Academic_Life(201808-Present)\2_Main_Research\CC_endo_regime\Program\Main_FCRET\Computing_code_uploaded\202510_NS_FCRGT') % Notebook 
% cd('C:\Users\user\Dropbox\Academic_Life(201808-Present)\2_Main_Research\CC_endo_regime\Program\Main_FCRET') % Office Desktop
%% Set-up
K_ind = 4; Non_dimX = (1:2); Non_dimY = (1:2);
%% Data Loading 
load('tic_g.mat','tic'); load('TPDN_g.mat','F'); GTemp_Grid = tic; GTemp0 = F(2:end,:); T = size(GTemp0,1);  
load('GRP_den2_full.mat','x_t','Fn'); GRP_Gr0 = Fn; GRP_Grid = x_t;  Date = (1951:1:2019); 
%% Data Visualization
%%% Figure 1
figure;
mesh(GTemp_Grid',(1951:2019),GTemp0)
xlabel('Land Temperature Anomaly','Fontsize',15,'fontweight','b'); 
ylabel('Year','Fontsize',35,'fontweight','b')
xlim([GTemp_Grid(1)-0.05 GTemp_Grid(end)+0.05])
ylim([1951 2019])
set(gca,'Linewidth',3.0,'box','off','Ticklength',[0 0],'Fontsize',35,'fontweight','b');
figure;
mesh(GRP_Grid,(1951:2019),GRP_Gr0)
xlabel('GRP Growth Rate','Fontsize',15,'fontweight','b'); 
ylabel('Year','Fontsize',35,'fontweight','b')
xlim([x_t(1)-0.05 x_t(end)+0.05])
ylim([1951 2019])
set(gca,'Linewidth',3.0,'box','off','Ticklength',[0 0],'Fontsize',35,'fontweight','b');
Sum_stat = KS_Desc_stat(GTemp_Grid',GTemp0); 
figure;
plot(Date,Sum_stat(:,1),'r','LineWidth',2.5); hold on;
plot(Date,sqrt(Sum_stat(:,2)),'b','LineWidth',2.5); hold on;
axis tight; grid on;
legend('Mean','Standard Deviation','Location','Northwest','Orientation','horizontal'); legend boxoff;
set(gca,'Linewidth',3.0,'box','on','Ticklength',[0 0],'Fontsize',35,'fontweight','b');
Sum_stat2 = KS_Desc_stat(GRP_Grid,GRP_Gr0); 
figure;
plot(Date,Sum_stat2(:,1),'r','LineWidth',2.5); hold on;
plot(Date,sqrt(Sum_stat2(:,2)),'b','LineWidth',2.5); hold on;
axis tight; grid on;
legend('Mean','Standard Deviation','Location','Northwest','Orientation','horizontal'); legend boxoff;
set(gca,'Linewidth',3.0,'box','on','Ticklength',[0 0],'Fontsize',35,'fontweight','b');
%% Centered-Log Ratio (CLR) Transformation & Nonstationary Decomposition
GTemp = zeros(size(GTemp0,1),size(GTemp0,2)); GRP_Gr = zeros(size(GRP_Gr0,1),size(GRP_Gr0,2)); 
for i=1:T
   GTemp(i,:) = log(GTemp0(i,:)/geomean(GTemp0(i,:)));  GRP_Gr(i,:) = log(GRP_Gr0(i,:)/geomean(GRP_Gr0(i,:))); 
end
[ N_GTemp, Inv_N_GTemp, S_GTemp1, S_GTemp, S_GTemp3, Inv_S_GTemp] = density_decom(GTemp_Grid, GTemp , max(Non_dimX));
[N_GRP_Gr,Inv_N_GRP_Gr,S_GRP_Gr1,S_GRP_Gr,S_GRP_Gr3,Inv_S_GRP_Gr] = density_decom( GRP_Grid', GRP_Gr, max(Non_dimY));
%% Temporal Demeaning Procedure
E_GTemp = mean(GTemp)';    GTemp = GTemp-kron(ones(T,1),E_GTemp'); 
E_GRP_Gr = mean(GRP_Gr)'; GRP_Gr = GRP_Gr-kron(ones(T,1),E_GRP_Gr'); 
%% Generating the Proposed Estimator
%%% Figure 2
Xraw_tem = GTemp0'; T1 = size(Xraw_tem, 2); indx = floor(T1 / 2); GW2_size = 5;

figure;
plot(GTemp_Grid,mean(Xraw_tem(:,1:indx),2),'k--','LineWidth',3.5); hold on; grid on;
plot(GTemp_Grid,mean(Xraw_tem(:,(indx+1):T1),2),'k','LineWidth',3.5); hold on;
axis([GTemp_Grid(1) GTemp_Grid(end) -0.01 0.45]); legend('1951-1984 Avg', '1985-2019 Avg', 'Location', 'Northwest','box','off');
title('Global Warming Function (GW1)','fontsize',35,'fontweight','b','LineWidth',10.0); 
xlabel('Temperature Anomaly','fontsize',35,'fontweight','b'); set(gca,'Linewidth',3.0,'box','on','Ticklength',[0 0],'Fontsize',35,'fontweight','b');
figure;
plot(GTemp_Grid,mean(Xraw_tem(:,1:GW2_size),2),'k--','LineWidth',3.5); hold on; grid on; 
plot(GTemp_Grid,mean(Xraw_tem(:,(T1-GW2_size+1):T1),2),'k','LineWidth',3.5); hold on;
axis([GTemp_Grid(1) GTemp_Grid(end) -0.01 0.45]); legend('1951-1955 Avg', '2015-2019 Avg', 'Location', 'Northwest','box','off');
title('Global Warming Function (GW2)','fontsize',35,'fontweight','b','LineWidth',10.0); 
xlabel('Temperature Anomaly','fontsize',35,'fontweight','b'); set(gca,'Linewidth',3.0,'box','on','Ticklength',[0 0],'Fontsize',35,'fontweight','b');

% Total-run Estimator
Tef0_gw0 = f_kappa_est_ci(0, GTemp_Grid, GTemp0, GTemp, GRP_Grid, GRP_Gr0, GRP_Gr, Inv_N_GTemp, Inv_S_GTemp, K_ind, Non_dimY, Non_dimX, 1);
Tef0_gw1 = f_kappa_est_ci(1, GTemp_Grid, GTemp0, GTemp, GRP_Grid, GRP_Gr0, GRP_Gr, Inv_N_GTemp, Inv_S_GTemp, K_ind, Non_dimY, Non_dimX, 1);
Tef0_gw0_GW2 = f_kappa_est_ci_GW(0, GTemp_Grid, GTemp0, GTemp, GRP_Grid, GRP_Gr0, GRP_Gr, Inv_N_GTemp, Inv_S_GTemp, K_ind, Non_dimY, Non_dimX, 1);
Tef0_gw1_GW2 = f_kappa_est_ci_GW(1, GTemp_Grid, GTemp0, GTemp, GRP_Grid, GRP_Gr0, GRP_Gr, Inv_N_GTemp, Inv_S_GTemp, K_ind, Non_dimY, Non_dimX, 1);
% Short-run Estimator
Sef0_gw0 = f_kappa_est_ci(0, GTemp_Grid, GTemp0, GTemp, GRP_Grid, GRP_Gr0, GRP_Gr, Inv_N_GTemp, Inv_S_GTemp, K_ind, Non_dimY, Non_dimX, 0);
Sef0_gw1 = f_kappa_est_ci(1, GTemp_Grid, GTemp0, GTemp, GRP_Grid, GRP_Gr0, GRP_Gr, Inv_N_GTemp, Inv_S_GTemp, K_ind, Non_dimY, Non_dimX, 0);
Sef0_gw0_GW2 = f_kappa_est_ci_GW(0, GTemp_Grid, GTemp0, GTemp, GRP_Grid, GRP_Gr0, GRP_Gr, Inv_N_GTemp, Inv_S_GTemp, K_ind, Non_dimY, Non_dimX, 0);
Sef0_gw1_GW2 = f_kappa_est_ci_GW(1, GTemp_Grid, GTemp0, GTemp, GRP_Grid, GRP_Gr0, GRP_Gr, Inv_N_GTemp, Inv_S_GTemp, K_ind, Non_dimY, Non_dimX, 0);

%%% Figure 3
index = 2:(length(GRP_Grid) - 1);
figure('WindowState','maximized'); 
subplot(1,2,1)
plot(GRP_Grid(index), Tef0_gw0(index,1), 'r','LineWidth',2.0); hold on;
plot(GRP_Grid(index), Tef0_gw1(index,1), 'k','LineWidth',4.0); hold on;
plot(GRP_Grid(index), Tef0_gw1(index,2), 'k--','LineWidth',3.0); hold on;
plot(GRP_Grid(index), Tef0_gw1(index,3), 'k--','LineWidth',3.0); hold on;
axis([min(GRP_Grid) max(GRP_Grid) -0.45 0.45]); grid on;
xlabel('GRP Growth Rate','fontsize',35,'fontweight','b');
title('Total-run Response Function','fontsize',30,'fontweight','b','LineWidth',10.0); 
legend('GW1 at \kappa=0','GW1 at \kappa=1', '95% CI', '','NumColumns', 1, 'Location', 'Southwest','box','off');
set(gca,'Linewidth',3.0,'box','on','Ticklength',[0 0],'Fontsize',30,'fontweight','b');
subplot(1,2,2) 
plot(GRP_Grid(index), Sef0_gw0(index,1), 'r','LineWidth',2.0); hold on;
plot(GRP_Grid(index), Sef0_gw1(index,1), 'k','LineWidth',4.0); hold on;
plot(GRP_Grid(index), Sef0_gw1(index,2), 'k--','LineWidth',3.0); hold on;
plot(GRP_Grid(index), Sef0_gw1(index,3), 'k--','LineWidth',3.0); hold on;
axis([min(GRP_Grid) max(GRP_Grid) -0.05 0.05]); grid on;
xlabel('GRP Growth Rate','fontsize',35,'fontweight','b');
title('Short-run Response Function','fontsize',30,'fontweight','b','LineWidth',10.0); 
legend('GW1 at \kappa=0','GW1 at \kappa=1', '95% CI', '','NumColumns', 1, 'Location', 'Southeast','box','off');
set(gca,'Linewidth',3.0,'box','on','Ticklength',[0 0],'Fontsize',30,'fontweight','b');
figure('WindowState','maximized'); 
subplot(1,2,1)
plot(GRP_Grid(index), Tef0_gw0_GW2(index,1), 'r','LineWidth',2.0); hold on;
plot(GRP_Grid(index), Tef0_gw1_GW2(index,1), 'k','LineWidth',4.0); hold on;
plot(GRP_Grid(index), Tef0_gw1_GW2(index,2), 'k--','LineWidth',3.0); hold on;
plot(GRP_Grid(index), Tef0_gw1_GW2(index,3), 'k--','LineWidth',3.0); hold on;
axis([min(GRP_Grid) max(GRP_Grid) -0.8 0.8]); grid on;
xlabel('GRP Growth Rate','fontsize',30,'fontweight','b');
title('Total-run Response Function','fontsize',35,'fontweight','b','LineWidth',10.0); 
legend('GW2 at \kappa=0','GW2 at \kappa=1', '95% CI', '','NumColumns', 1, 'Location', 'Southwest','box','off');
set(gca,'Linewidth',3.0,'box','on','Ticklength',[0 0],'Fontsize',30,'fontweight','b');
subplot(1,2,2) 
plot(GRP_Grid(index), Sef0_gw0_GW2(index,1), 'r','LineWidth',2.0); hold on;
plot(GRP_Grid(index), Sef0_gw1_GW2(index,1), 'k','LineWidth',4.0); hold on;
plot(GRP_Grid(index), Sef0_gw1_GW2(index,2), 'k--','LineWidth',3.0); hold on;
plot(GRP_Grid(index), Sef0_gw1_GW2(index,3), 'k--','LineWidth',3.0); hold on;
axis([min(GRP_Grid) max(GRP_Grid) -0.15 0.2]); grid on;
xlabel('GRP Growth Rate','fontsize',30,'fontweight','b');
title('Short-run Response Function','fontsize',35,'fontweight','b','LineWidth',10.0); 
legend('GW2 at \kappa=0','GW2 at \kappa=1', '95% CI', '','NumColumns', 1, 'Location', 'Southwest','box','off');
set(gca,'Linewidth',3.0,'box','on','Ticklength',[0 0],'Fontsize',30,'fontweight','b');

%% Climate Change Impacts Projection
Wght = linspace(0,1.5,30); Wght(7) = 0.3; Wght(13) = 0.6; Wght(19) = 0.9; Wght(24) = 1.2; Wght(30) = 1.5; 
[Target0_grid0, ref_recover, Proj_raw1, Proj_raw2, Proj_raw3, Proj_raw4,Proj_raw5,Sum_stat] = CC_Frac_Impacts(GRP_Grid,GRP_Gr0,Tef0_gw1,Wght);
[Target0_grid0_GW2, ref_recover_GW2, Proj_raw1_GW2, Proj_raw2_GW2, Proj_raw3_GW2, Proj_raw4_GW2,Proj_raw5_GW2,Sum_stat_GW2] = CC_Frac_Impacts(GRP_Grid,GRP_Gr0,Tef0_gw1_GW2,Wght);
%%% Figure 4
figure('WindowState','maximized');
plot(Target0_grid0(2:end-1), ref_recover(2:end-1), 'k-.', 'LineWidth', 4.5); hold on;
plot(Target0_grid0(2:end-1), Proj_raw1(2:end-1), 'k:', 'LineWidth', 3.0); hold on;
plot(Target0_grid0(2:end-1), Proj_raw2(2:end-1), 'k:', 'LineWidth', 3.0); hold on;
plot(Target0_grid0(2:end-1), Proj_raw3(2:end-1), 'k:', 'LineWidth', 3.0); hold on;
plot(Target0_grid0(2:end-1), Proj_raw4(2:end-1), 'k:', 'LineWidth', 3.0); hold on;
plot(Target0_grid0(2:end-1), Proj_raw5(2:end-1), 'k', 'LineWidth', 4.5); hold on;
xlabel('GRP Growth Rate','fontsize',35,'fontweight','b'); 
axis([min(Target0_grid0) max(Target0_grid0) 0.2 1.7]); grid on;
legend('Avg GRP(51-84)','q=0.3', 'q=0.6', 'q=0.9', 'q=1.2','q=1.5', 'Location', 'Northoutside','box','off','Orientation','horizontal','Fontsize',40,'fontweight','b');
set(gca,'Linewidth',3.0,'box','on','Ticklength',[0 0],'Fontsize',35,'fontweight','b');

figure('WindowState','maximized');
plot(Target0_grid0(2:end-1), ref_recover_GW2(2:end-1), 'k-.', 'LineWidth', 4.5); hold on;
plot(Target0_grid0(2:end-1), Proj_raw1_GW2(2:end-1), 'k:', 'LineWidth', 3.0); hold on;
plot(Target0_grid0(2:end-1), Proj_raw2_GW2(2:end-1), 'k:', 'LineWidth', 3.0); hold on;
plot(Target0_grid0(2:end-1), Proj_raw3_GW2(2:end-1), 'k:', 'LineWidth', 3.0); hold on;
plot(Target0_grid0(2:end-1), Proj_raw4_GW2(2:end-1), 'k:', 'LineWidth', 3.0); hold on;
plot(Target0_grid0(2:end-1), Proj_raw5_GW2(2:end-1), 'k', 'LineWidth', 4.5); hold on;
xlabel('GRP Growth Rate','fontsize',35,'fontweight','b'); 
axis([min(Target0_grid0) max(Target0_grid0) 0.2 1.7]); grid on;
legend('Avg GRP(51-84)','q=0.3', 'q=0.6', 'q=0.9', 'q=1.2','q=1.5', 'Location', 'Northoutside','box','off','Orientation','horizontal','Fontsize',40,'fontweight','b');
set(gca,'Linewidth',3.0,'box','on','Ticklength',[0 0],'Fontsize',35,'fontweight','b');

figure('WindowState','maximized');
subplot(2,1,1)
plot(Wght, Sum_stat(2:end,1), 'b', 'LineWidth', 2.5); hold on;
plot(Wght, Sum_stat_GW2(2:end,1), 'r', 'LineWidth', 2.5); hold on;
title('Mean','fontsize',35,'fontweight','b','LineWidth',10.0); axis tight; grid on;
legend('GW1', 'GW2', 'Location', 'Northeast','box','off','Orientation','horizontal');
set(gca,'Linewidth',3.0,'box','on','Ticklength',[0 0],'Fontsize',33,'fontweight','b');
subplot(2,1,2)
plot(Wght, Sum_stat(2:end,2), 'b', 'LineWidth', 2.5); hold on;
plot(Wght, Sum_stat_GW2(2:end,2), 'r', 'LineWidth', 2.5); hold on;
legend('GW1', 'GW2', 'Location', 'Northwest','box','off','Orientation','horizontal');
title('Variance','fontsize',35,'fontweight','b','LineWidth',10.0); axis tight; grid on;
set(gca,'Linewidth',3.0,'box','on','Ticklength',[0 0],'Fontsize',33,'fontweight','b');

% EOP

