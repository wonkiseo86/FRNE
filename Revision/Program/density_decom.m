function [N_Fn,Inv_N_Fn, Est_S_Fg1, Est_S_Fg2, Est_S_Fg3,Inv_Est_S_Fg2] = density_decom(X_grid, X, Non_dim)
% X: CLR transformed density, Non_dim: non-stationary dimension of X

Figure_on = 1;

T = size(X,1);
Qhat = X'*X; % Estimated Covariance Operator
[eigvecs, eigvals] = eig(Qhat); % Eigenvalues and Eigenvectors
 [eigvals, indsrt] = sort(diag(eigvals),'descend'); % Ordered eigenvalues
eigvecs = eigvecs(:,indsrt); % Ordered eigenvectors

N_Fn = zeros(size(X,1),size(X,2)); S_Fn = zeros(size(X,1),size(X,2));
  X1 = eigvecs(:,1:Non_dim);         X2 = eigvecs(:,Non_dim+1:end);
PI_N = X1*(X1'*X1)^(-1)*X1';       PI_S = X2*(X2'*X2)^(-1)*X2'; 

for i = 1:T
    N_Fn(i,:)=PI_N*X(i,:)'; % Nonstationary component 
    S_Fn(i,:)=PI_S*X(i,:)'; % Stationary component 
end

Est_S_Fg1 = X - N_Fn;  % Get stationary component by subtracting nonstationary component from raw density
Est_S_Fg2 = S_Fn; % Get stationary component with stationary component

% for i = 1:T
%     temp0 = zeros(length(eigvecs),1);
%     for j=1:length(eigvecs)-2
%     temp = (X2(:,j)'*X(i,:)')*X2(:,j); % Stationary component
%     temp0 = temp0+temp;
%     end
%     S_Fn(i,:) = temp0;
% end

Est_S_Fg3 = S_Fn;


if Figure_on == 1
    x_t = X_grid;
    %%% Inverse CLR transformation
    Inv_N_Fn = zeros(size(N_Fn,1),size(N_Fn,2)); Inv_Est_S_Fg2 = zeros(size(Est_S_Fg2,1),size(Est_S_Fg2,2));
    for i =1:size(N_Fn,1)
    Temp1 = mrsum(x_t',exp(N_Fn(i,:))'); Temp2 = mrsum(x_t',exp(Est_S_Fg2(i,:))');
    Inv_N_Fn(i,:) = exp(N_Fn(i,:))/Temp1; Inv_Est_S_Fg2(i,:) = exp(Est_S_Fg2(i,:))/Temp2;
    end
% figure;
% subplot(1,2,1)
% mesh(x_t,(1951:1951+T-1),Inv_N_Fn)
% % xlabel('Temperature Anomaly','Fontsize',25,'fontweight','b'); 
% xlabel('GRP Growth Rate','Fontsize',15,'fontweight','b'); 
% title('Nonstationary Component','fontsize',25,'fontweight','b')
% ylabel('Year','Fontsize',25,'fontweight','b')
% xlim([x_t(1)-0.05 x_t(end)+0.05])
% ylim([1951 1951+T-1])
% set(gca,'Linewidth',3.0,'box','off','Ticklength',[0 0],'Fontsize',25,'fontweight','b');
% subplot(1,2,2)
% mesh(x_t,(1951:1951+T-1),Inv_Est_S_Fg2)
% % xlabel('Temperature Anomaly','Fontsize',25,'fontweight','b'); 
% xlabel('GRP Growth Rate','Fontsize',15,'fontweight','b'); 
% title('Stationary Component','fontsize',25,'fontweight','b')
% ylabel('Year','Fontsize',25,'fontweight','b')
% xlim([x_t(1)-0.05 x_t(end)+0.05])
% ylim([1951 1951+T-1])
% set(gca,'Linewidth',3.0,'box','off','Ticklength',[0 0],'Fontsize',25,'fontweight','b');
% % grid on;
    
s=x_t'; sumstats = [];  sumstats2 = []; sy=1951; ey=1951+T-1;

for i=1:T
    f = Inv_N_Fn(i,:)';
    ff= s.*f ; mm= mrsum(s,ff); ss=s;
    for k=2:4; ss = [ss (s-mm).^k]; end
    ff2 = ss.*repmat(f,1,size(ss,2)) ; a0  = [(sy+i-1) mrsum(s,ff2)'] ; sumstats = [sumstats; a0];
end
N_c1 = sumstats(:,2); N_c2 = sumstats(:,3); N_c3 = sumstats(:,4)./(N_c2.^(3/2)) ; N_c4 = sumstats(:,5)./(N_c2.^2) ; 

for i=1:T
    f = Inv_Est_S_Fg2(i,:)';
    ff= s.*f ; mm= mrsum(s,ff); ss=s;
    for k=2:4; ss = [ss (s-mm).^k]; end
    ff2 = ss.*repmat(f,1,size(ss,2)) ; a0  = mrsum(s,ff2)'; sumstats2 = [sumstats2; a0];
end
S_c1 = sumstats2(:,1); S_c2 = sumstats2(:,2); S_c3 = sumstats2(:,3)./(S_c2.^(3/2)) ; S_c4 = sumstats2(:,4)./(S_c2.^2) ; 

% figure;
% subplot(2,2,1)
% plot((sy:ey),N_c1,'-r','LineWidth',2.0); hold on;
% plot((sy:ey),S_c1,'-b','LineWidth',2.5); hold on;
% xlim([sy ey])
% axis tight;
% title('Mean','fontsize',25,'fontweight','b')
% set(gca,'Linewidth',3.0,'box','off','Ticklength',[0 0],'Fontsize',20,'fontweight','b');
% 
% subplot(2,2,2)
% plot((sy:ey),N_c2,'-r','LineWidth',2.0); hold on;
% plot((sy:ey),S_c2,'-b','LineWidth',2.5); hold on;
% xlim([sy ey])
% axis tight;
% title('Variance','fontsize',25,'fontweight','b')
% set(gca,'Linewidth',3.0,'box','off','Ticklength',[0 0],'Fontsize',20,'fontweight','b');
% 
% subplot(2,2,3)
% plot((sy:ey),N_c3,'-r','LineWidth',2.0); hold on;
% plot((sy:ey),S_c3,'-b','LineWidth',2.5); hold on;
% xlim([sy ey])
% axis tight;
% xlabel('off','color','white','LineWidth',2.0)
% title('Skewness','fontsize',25,'fontweight','b')
% legend('Nonstationary Component','Stationary Component','Location','SouthOutside','Orientation','horizontal'); legend boxoff;
% set(gca,'Linewidth',3.0,'box','off','Ticklength',[0 0],'Fontsize',20,'fontweight','b');
% 
% subplot(2,2,4)
% plot((sy:ey),N_c4,'-r','LineWidth',2.0); hold on;
% plot((sy:ey),S_c4,'-b','LineWidth',2.5); hold on; 
% axis tight;
% xlim([sy ey])
% xlabel('off','color','white','LineWidth',2.0)
% title('Kurtosis','fontsize',25,'fontweight','b')
% legend('Nonstationary Component','Stationary Component','Location','SouthOutside','Orientation','horizontal'); legend boxoff;
% set(gca,'Linewidth',3.0,'box','off','Ticklength',[0 0],'Fontsize',20,'fontweight','b');

end





