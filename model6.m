
% ~~~~code for fitting UPS data~~~~~~~~ 

data_1 = load("034_1.0331_35MLAg.txt"); % load data file1..substrate
data_2 = load("034_1.0331_HATCN.txt"); % load data file2..adsorbed system
l1 = length(data_1(:,1));
l2 = length(data_2(:,1));
EF = 22.08; %Fermi level KE

%substrate
x1 = data_1(:,1);
y1 = data_1(:,2);
% x1(:,1) = x1(:,1)-EF; %KE to BE

%adsorbed system
x2 = data_2(:,1);
y2 = data_2(:,2);
% plot(x1,y1,x2,y2)
% x2(:,1) = x2(:,1)-EF; %KE to BE
BE_start=125;

x2_ad = x2(BE_start:end,1);
y2_ad = y2(BE_start:end,1);
y2_ad = y2_ad-min(y2_ad(:,1));
% y2_ad =smoothdata(y2_ad);
%substrate adjusted for adsorbed data
x1_ad = x1(BE_start:end,1);
y1_ad = y1(BE_start:end,1);
y1_ad = y1_ad-min(y1_ad(:,1));
% y1_ad =smoothdata(y1_ad);
% %changing the EF energy
x1_ad=x1_ad-0.04;
x2_ad=x2_ad-0.04;

plot(x1_ad,y1_ad,x2_ad,y2_ad)
xlabel('BE (eV)');
ylabel('counts');
legend('Ag','1 ML HATCN','Location','northeast')



% ~~Tougaard BG~~~
%~~~~Tougaard BG for substrate
Tgd_sub = zeros(length(x1_ad(:,1)),1);
dE = abs(x1_ad(1,1)-x1_ad(2,1));
P = 20;
for p = 1:length(x1_ad(:,1))
    for q = p:length(x1_ad(:,1))
        Tgd_sub(p,1) = Tgd_sub(p,1) + y1_ad(q,1)*dE*((x1_ad(q,1)-x1_ad(p,1))/(P+(x1_ad(q,1)-x1_ad(p,1).^2).^2));
    end
end
B1 = y1_ad(1,1)/Tgd_sub(1,1);
Tgd_sub = B1*Tgd_sub;
plot(x1_ad,Tgd_sub,x1_ad,y1_ad)
% plot(x1_ad,y1_ad-Tgd_sub)
% hold on
%~~~~Tougaard BG for adsorbed system
 Tgd_ads = zeros(length(x2_ad(:,1)),1);
dE = abs(x2_ad(1,1)-x2_ad(2,1));
P =20;
% P = 50;
for p = 1:length(x2_ad(:,1))
    for q = p:length(x2_ad(:,1))
        Tgd_ads(p,1) = Tgd_ads(p,1) + y2_ad(q,1)*dE*((x2_ad(q,1)-x2_ad(p,1))/(P+(x2_ad(q,1)-x2_ad(p,1).^2).^2));
    end
end
B2 = y2_ad(1,1)/Tgd_ads(1,1);
% plot(x2_ad,B2*Tgd_ads,x2_ad,y2_ad)
% hold on


% plot(x1_ad,y1_ad,x1_ad,y1_ad-B1*Tgd_sub)
% plot(x2_ad,y2_ad-B2*Tgd_ads)
% legend('Ag','HATCN/Ag')
% hold off

%finding efffective temp from Ag spectrum
start_param = [400,5,-300]; 
lower_param = [0,0,-1E5];
upper_param = [1E5,150,1E5];

ftopt = fitoptions('Method','NonlinearLeastSquares','StartPoint',start_param,'Lower',lower_param,'Upper',upper_param);
ftmethod = fittype('K*((1+ exp(x1_ad.*1000/(25+M))).^(-1))+F*x1_ad','independent',{'x1_ad'},'dependent',{'y1_ad'},'coefficients',{'K','M','F'},'options',ftopt);
% ftmethod = fittype('A*exp(-(((x1_ad-b).^2)/(c.^2)))+ D*(y1_ad) + E*Tgd_ads','problem',{'y1_ad','Tgd_ads'},'independent',{'x1_ad'},'dependent',{'y2'},'coefficients',{'A','b','c','D','E'},'options',ftopt);

[f,gof] = fit(x1_ad,y1_ad,ftmethod);
cf = coeffvalues(f);
r2 = gof.rsquare;

FD = zeros();
lnBG = zeros();
f = zeros();
for i=1:length(x1_ad(:,1))
    FD(i,1) = cf(1)*((1+ exp(x1_ad(i,1).*1000/(25+cf(2)))).^(-1));
    lnBG(i,1)= cf(3)*x1_ad(i,1);

    f(i,1) = FD(i,1)+lnBG(i,1);

end
plot(x1_ad,y1_ad,x1_ad,f,'O',x1_ad,FD,x1_ad,lnBG)
T_eff = cf(2)+25;

%fitting HATCN spectrum

% start_param = [max(y2_ad(:,1)),-0.05,0.2,0]; 
% lower_param = [0,-0.5,0,0];
% upper_param = [4*max(y2_ad(:,1)),0.2,0.4,1];
C=0.35;

start_param = [max(y2_ad(:,1)),0.05,0.02,0.3]; 
lower_param = [0,-0.5,-C,0];
upper_param = [2*max(y2_ad(:,1)),0.25,0.15,1];
T_eff =40.07;%manually choose Fermi broadening

ftopt = fitoptions('Method','NonlinearLeastSquares','StartPoint',start_param,'Lower',lower_param,'Upper',upper_param);
% ftmethod = fittype('A*exp(-(((x2_ad-b).^2)/(c.^2)))*((1+ exp(x2_ad.*1000/(T_eff))).^(-1))+ D*(y1_ad-Tgd_sub) + B2*Tgd_ads','problem',{'T_eff','y1_ad','Tgd_sub','B2','Tgd_ads'},'independent',{'x2_ad'},'dependent',{'y_fit'},'coefficients',{'A','b','c','D'},'options',ftopt);
ftmethod = fittype('A*exp(-(((x2_ad-b).^2)/((C+e).^2)))*((1+ exp(x2_ad.*1000/(T_eff))).^(-1))+ D*(y1_ad-Tgd_sub) + B2*Tgd_ads','problem',{'C','T_eff','y1_ad','Tgd_sub','B2','Tgd_ads'},'independent',{'x2_ad'},'dependent',{'y2_ad'},'coefficients',{'A','b','e','D'},'options',ftopt);


% [f,gof,output] = fit(x2_ad,y_fit,ftmethod,'problem',{T_eff,y1_ad,Tgd_sub,B2,Tgd_ads});
[f,gof,output] = fit(x2_ad,y2_ad,ftmethod,'problem',{C,T_eff,y1_ad,Tgd_sub,B2,Tgd_ads});

cf = coeffvalues(f);
pk = cf(2);
HW = C+cf(3);
r2 = gof.rsquare;

%~~~plot~~~

fitfn = zeros();
LUMO = zeros();
sp_band = zeros();
diff_y = zeros();
Ag_BG = zeros(); 
lnBG = zeros();
Tgd_BG = zeros();
BE_final = 1;
step = 0.01;
x2_ad_extended = transpose(linspace(x2_ad(1,1),BE_final,(BE_final-x2_ad(1,1))/step));
for i=1:length(x2_ad_extended(:,1))
    LUMO(i,1) = cf(1)*exp(-(((x2_ad_extended(i,1)-cf(2)).^2)/((C+cf(3)).^2)));

end
for i=1:length(x2_ad(:,1))
%     LUMO(i,1) = cf(1)*exp(-(((x2_ad(i,1)-cf(2)).^2)/(cf(3).^2)));
    LUMO(i,1) = cf(1)*exp(-(((x2_ad(i,1)-cf(2)).^2)/((C+cf(3)).^2)));

    FD(i,1) = ((1+ exp(x2_ad(i,1).*1000/(T_eff))).^(-1));
    Ag_BG(i,1)= cf(4)*(y1_ad(i,1)-Tgd_sub(i,1));
%     lnBG(i,1)= (cf(5)*x2_ad(i,1));
    Tgd_BG(i,1)= (B2*Tgd_ads(i,1));
    fitfn(i,1) = LUMO(i,1)*FD(i,1)+ Ag_BG(i,1)+Tgd_BG(i,1);

end

area(x2_ad_extended,LUMO,'FaceColor','c','FaceAlpha',.3,'EdgeAlpha',.1)
xlim([x2_ad_extended(1,1) x2_ad_extended(end,1)])
ylim([0 max(y2_ad)])
hold on
plot(x2_ad,fitfn,'or',x2_ad,y2_ad,x2_ad,Ag_BG,x2_ad,Tgd_BG,'MarkerSize',3,'LineWidth',1)
legend('G(E)','fit','HATCN/35 ML Ag/Cu(111)','Ag sp','Tougaard BG','Location','west')
xlabel('Binding Energy(eV)')
ylabel('Counts(a.u.)')
a = get(gca,'XTickLabel'); 
set(gca,'fontsize',18,'FontWeight','bold');
b = get(gca,'XTick');
set(gca,'LineWidth',2);
c = get(gca,'legend'); 
set(gca,'fontsize',14,'FontWeight','bold');
legend('boxoff')
hold off
% box on
% plot(x2_ad,y2_ad,'ok','MarkerSize',6)
% hold on
% % plot(x2_ad_extended,LUMO,'or','MarkerSize',3)
% % hold on
% area(x2_ad_extended,LUMO,'FaceColor','c','FaceAlpha',.3,'EdgeAlpha',.1)
% xlim([-1 1])
% % legend('HATCN/1ML Ag/Cu(111)','fit','G(E)','Location','southeast')
% xlabel('Binding Energy(eV)')
% set(gca,'ytick',[])
% xline(pk,'--r','LineWidth',3)
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
