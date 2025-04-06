data = load('FM spin down inv.txt');
V = data(:,1);
DOS = data(:,2);

% c = min(DOS(:,1));
% DOS(:,1) = DOS(:,1)-c;
x = V(:,1);
y = DOS(:,1);
plot(x,y);

%initial guess fit
[param,~] = fit(x,y,'gauss2');
pval = coeffvalues(param);
plot(param,x,y);
m = ((y(end,1)-y(1,1))/(x(end,1)-x(1,1)));
K = y(1,1);
% Fit w/ 2 Gaussians
% start_param = [pval(1,1),pval(1,2),pval(1,3),pval(1,5),pval(1,6),m,0];
% % start_param = [pval(1,1),0.09,0.078,0.026,0.1,m,0];
% lower_param = [0,V(1,1),0.05,V(1,1),0.05,0,0];
% % lower_param = [0,V(1,1),0.065,V(1,1),0.065,-Inf,0];
% upper_param = [+Inf,V(end,1),0.1,V(end,1),0.1,+Inf,max(y(:,1))];
% 
% 
% 
% ftopt = fitoptions('Method','NonlinearLeastSquares','Lower',lower_param,'Upper',upper_param,'StartPoint',start_param,'TolX',1e-6,'TolFun',1e-8,'MaxIter',800);
% ftmethod = fittype('a1*exp(-((x-b1)/c1)^2) + a1*exp(-((x-b2)/c2)^2) + D*x + E','independent',{'x'},'dependent',{'y'},'coefficients',{'a1','b1','c1','b2','c2','D','E'},'options',ftopt);

start_param = [pval(1,1),pval(1,2),pval(1,3),pval(1,5),m,0];
lower_param = [0,V(1,1),0,V(1,1),0,0];
upper_param = [+Inf,V(end,1),0.10,V(end,1),+Inf,max(y(:,1))];

x = V(:,1);
y = DOS(:,1);

ftopt = fitoptions('Method','NonlinearLeastSquares','Lower',lower_param,'Upper',upper_param,'StartPoint',start_param,'TolX',1e-6,'TolFun',1e-6,'MaxIter',800);
ftmethod = fittype('a1*exp(-((x-b1)/c1)^2) + a1*exp(-((x-b2)/c1)^2) + D*x + E','independent',{'x'},'dependent',{'y'},'coefficients',{'a1','b1','c1','b2','D','E'},'options',ftopt);

[f,gof,output] = fit(x,y,ftmethod);
cf = coeffvalues(f);
r2 = gof.rsquare;
figure
plot(f,x,y)
xlabel('Bias(V)','fontweight','bold','fontsize',22);
ylabel('PDOS(a.u.)','fontweight','bold','fontsize',22);
% xlim([-0.04 0.08]);
% xlim([-0.05 0.03]);
% ylim([0.7 1.2]);
a = get(gca,'XTickLabel'); 
set(gca,'fontsize',18,'FontWeight','bold');
b = get(gca,'XTick');
set(gca,'LineWidth',2);
legend('HATCN PDOS','fit','Location','northwest')
legend('boxoff')
box on

a1 = cf(1,1);
b1 = cf(1,2);
c1 = cf(1,3);
b2 = cf(1,4);
% c2 = cf(1,5);
% D = cf(1,6);
% E = cf(1,7);
D = cf(1,5);
E = cf(1,6);
t = zeros();
f1 = zeros();
f2 = zeros();
f3 = zeros();
for i=1:length(x(:,1))
    f1(i,1) = a1*exp(-((x(i,1)-b1)/c1)^2);
    f2(i,1) = a1*exp(-((x(i,1)-b2)/c1)^2);
    f3(i,1) = D*x(i,1) + E;
    t(i,1) = f1(i,1)+f2(i,1)+f3(i,1);
end
% T = table(V,t);
figure
% plot(x,y,'-o',V,t,V,f1,V,f2,'Linewidth',1)
plot(x,y,'-o',x,t,x,f1,x,f2,x,f3,'Linewidth',1)
xlabel('Bias(V)','fontweight','bold','fontsize',22);
ylabel('PDOS(a.u.)','fontweight','bold','fontsize',22);
xlim([x(1,1) x(end,1)]);
% xlim([-0.05 0.03]);
% ylim([0.7 1.2]);
a = get(gca,'XTickLabel'); 
set(gca,'fontsize',18,'FontWeight','bold');
b = get(gca,'XTick');
set(gca,'LineWidth',2);
legend('HATCN PDOS','fit','Gauss1','Gauss2','Location','northwest')
legend('boxoff')
box on

%Covariance matrix
J = output.Jacobian;
R = output.residuals;
N = output.numobs;
p = output.numparam;
MSE = (R'*R)/(N-p);
VC = inv(J'*J)*MSE; %Covariance matrix

% Correlation coeficients
Corr = zeros();
for i = 1:p
    for j = 1:p

        Corr(i,j) = VC(i,j)/sqrt(VC(i,i)*VC(j,j));
    end
end
