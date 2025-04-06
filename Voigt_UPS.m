%Voigt fitting of UPS data with same center for  Gaussian and Lorentzian
%parts
data_1 = load("1ML_Ag"); % load data file1..substrate
data_2 = load("05_AR350degC2CAE100_22MLHCNAgCu111"); % load data file2..adsorbed system
l1 = length(data_1(:,1));
l2 = length(data_2(:,1));
EF = 22.08; %Fermi level KE

%substrate
x1 = data_1(:,1);
y1 = data_1(:,2);

%adsorbed system
x2 = data_2(:,1);
y2 = data_2(:,2);

BE_start=1;

x2_ad = x2(BE_start:end,1);
x2_ad = x2_ad-EF;
y2_ad = y2(BE_start:end,1);
y2_ad = y2_ad-min(y2_ad(:,1));

%substrate adjusted for adsorbed data
x1_ad = x1(BE_start:end,1);
x1_ad = x1_ad-EF;
y1_ad = y1(BE_start:end,1);
y1_ad = y1_ad-min(y1_ad(:,1));

%normalize both data to either max intensity or to BG
y1_ad = y1_ad/max(y1_ad(:,1));
y2_ad = y2_ad/max(y2_ad(:,1));

% figure
% plot(x1_ad,y1_ad,'-b',x2_ad,y2_ad,'-r','LineWidth',2)
% legend('22ML Ag','0.4 ML HATCN/22ML Ag','Location','northeast')
% xlabel('Binding Energy(eV)')
% ylabel('Counts(a.u.)')
% % xlim([-0.3 0.3])
% % ylim([0 1.2])
% % xline(pk1,'-r')
% a = get(gca,'XTickLabel'); 
% set(gca,'fontsize',18,'FontWeight','bold');
% b = get(gca,'XTick');
% set(gca,'LineWidth',2);
% c = get(gca,'legend');
% set(gca,'fontsize',14);
% legend('boxoff')
% box on

% ~~Tougaard BG~~~
%~~~~Tougaard BG substrate~~~~ 
Tgd_sub = zeros(length(x1_ad(:,1)),1);
dE = abs(x1_ad(1,1)-x1_ad(2,1));
P =21;
for p = 1:length(x1_ad(:,1))
    for q = p:length(x1_ad(:,1))
        Tgd_sub(p,1) = Tgd_sub(p,1) + y1_ad(q,1)*dE*((x1_ad(q,1)-x1_ad(p,1))/(P+(x1_ad(q,1)-x1_ad(p,1).^2).^2));
    end
end
B1 = y1_ad(1,1)/Tgd_sub(1,1);
Tgd_sub = B1*Tgd_sub;
plot(x1_ad,Tgd_sub,x1_ad,y1_ad)

%~~~~Tougaard BG adsorbed~~~~ 
Tgd_ads = zeros(length(x2_ad(:,1)),1);
dE = abs(x2_ad(1,1)-x2_ad(2,1));
P =21;
for p = 1:length(x2_ad(:,1))
    for q = p:length(x2_ad(:,1))
        Tgd_ads(p,1) = Tgd_ads(p,1) + y2_ad(q,1)*dE*((x2_ad(q,1)-x2_ad(p,1))/(P+(x2_ad(q,1)-x2_ad(p,1).^2).^2));
    end
end
B2 = y2_ad(1,1)/Tgd_ads(1,1);
Tgd_ads = B2*Tgd_ads;
% plot(x2_ad,Tgd_ads,x2_ad,y2_ad,x1_ad,Tgd_sub,x1_ad,y1_ad)
% plot(x2_ad,Tgd_ads,x2_ad,y2_ad)
[BE_Ag,loc1]= max(y1_ad(:,1));
Mx_Ag = x1_ad(loc1,1); %maxima position of substrate UPS
% [BE_HCN,loc2]= max(y2_ad(:,1));
% Mx_HCN = x2_ad(loc2,1); %maxima position of HATCN/Ag UPS
% pk_shift = Mx_HCN - Mx_Ag;
%~~~~~~~~~fitting using Voigt*Fermi-Dirac + Tougard BG~~~~~~~~~~~~~~
 
% substrate fitting
start_param = [1,-0.1,0.02,1,0.02,15,1]; 
lower_param = [0,-0.5,0,0,0,0,0];
upper_param = [1E5,0.25,0.09,1E3,1,150,1];

ftopt = fitoptions('Method','NonlinearLeastSquares','StartPoint',start_param,'Lower',lower_param,'Upper',upper_param);
% ftmethod = fittype('conv(F*exp(-(((x1_ad-G).^2)/((H).^2))),A*(C/((x1_ad-b).^2+((C).^2))),"same")*((1+ exp(x1_ad.*1000/(25+M))).^(-1))+ D*Tgd_sub','problem',{'Tgd_sub'},'independent',{'x1_ad'},'dependent',{'y1_ad'},'coefficients',{'F','G','H','A','C','b','M','D'},'options',ftopt);
%fit model with same center for both lorentz and gaussian
ftmethod = fittype('conv(F*exp(-(((x1_ad-G).^2)/((H).^2))),A*(C/((x1_ad-G).^2+((C).^2))),"same")*((1+ exp(x1_ad.*1000/(25+M))).^(-1))+ D*Tgd_sub','problem',{'Tgd_sub'},'independent',{'x1_ad'},'dependent',{'y1_ad'},'coefficients',{'F','G','H','A','C','M','D'},'options',ftopt);

[f,gof,output] = fit(x1_ad,y1_ad,ftmethod,'problem',{Tgd_sub});

cf = coeffvalues(f);
pk1 = cf(2);

HW1 = 0.5*cf(3);
HW2 = 0.5*cf(5);
r2 = gof.rsquare;
T_eff = 25+cf(6);

%plotting the output
fitfn = zeros();
fitfn2 = zeros();
G = zeros();
L= zeros();
FD = zeros();
BG = zeros();
vgt = zeros();
for i=1:length(x1_ad(:,1))
    G(i,1) = cf(1)*exp(-(((x1_ad(i,1)-cf(2)).^2)/((cf(3)).^2)));
    L(i,1) = cf(4)*(cf(5)/((x1_ad(i,1)-cf(2)).^2+((cf(5)).^2)));
    FD(i,1) = ((1+ exp(x1_ad(i,1).*1000/(25+cf(6)))).^(-1));
    BG(i,1) = cf(7)*Tgd_sub(i,1);
%     vgt(i,1) = conv(G(i,1),L(i,1),"same")*FD(i,1);

end
fitfn = conv(G,L,"same").*FD+BG;
% fitfn2 = vgt;
% q = cumtrapz(vgt);
FW_v = 0.5346*cf(5) +sqrt(0.2166*(cf(5)).^2 +(cf(3)).^2); %FWHM of voigt profile from literature
HW_v = 0.5*FW_v; %HWHM of voigt profile
C_v = cf(2); %center of Voigt profile from literature
[BE_v,loc]= max(fitfn(:,1));
Mx_v = x1_ad(loc,1); %maxima position of voigt function
param_Ag = [C_v HW_v T_eff]; %save Voigt parameters to an array

figure
plot(x1_ad,y1_ad,'ob',x1_ad,fitfn,'-r','LineWidth',2)
legend('10 ML Ag','fit','Location','northwest')
xlabel('Binding Energy(eV)')
ylabel('Counts(a.u.)')
% xlim([-0.3 0.3])
ylim([0 1.2])
% xline(pk1,'-r')
a = get(gca,'XTickLabel'); 
set(gca,'fontsize',22,'FontWeight','bold');
b = get(gca,'XTick');
set(gca,'LineWidth',2);
c = get(gca,'legend');
set(gca,'fontsize',16);
legend('boxoff')
box on

% %adsorbate fitting
% 
% start_param = [1,-0.025,0.02,1,0.009,5,1]; 
% lower_param = [0,-0.06,0,0,0,0,0];
% upper_param = [1E3,0,1,1E3,0.1,150,1];
% 
% ftopt = fitoptions('Method','NonlinearLeastSquares','StartPoint',start_param,'Lower',lower_param,'Upper',upper_param);
% % ftmethod = fittype('conv(F*exp(-(((x1_ad-G).^2)/((H).^2))),A*(C/((x1_ad-b).^2+((C).^2))),"same")*((1+ exp(x1_ad.*1000/(25+M))).^(-1))+ D*Tgd_sub','problem',{'Tgd_sub'},'independent',{'x1_ad'},'dependent',{'y1_ad'},'coefficients',{'F','G','H','A','C','b','M','D'},'options',ftopt);
% %fit model with same center for both lorentz and gaussian
% ftmethod = fittype('conv(F*exp(-(((x2_ad-G).^2)/((H).^2))),A*(C/((x2_ad-G).^2+((C).^2))),"same")*((1+ exp(x2_ad.*1000/(25+M))).^(-1))+ D*Tgd_ads','problem',{'Tgd_ads'},'independent',{'x2_ad'},'dependent',{'y2_ad'},'coefficients',{'F','G','H','A','C','M','D'},'options',ftopt);
% 
% [f,gof,output] = fit(x2_ad,y2_ad,ftmethod,'problem',{Tgd_ads});
% 
% cf = coeffvalues(f);
% pk1 = cf(2);
% 
% HW1 = 0.5*cf(3);
% HW2 = 0.5*cf(5);
% r2 = gof.rsquare;
% T_eff = 25+cf(6);
% 
% %plotting the output
% fitfn = zeros();
% fitfn2 = zeros();
% G = zeros();
% L= zeros();
% FD = zeros();
% BG = zeros();
% vgt = zeros();
% for i=1:length(x2_ad(:,1))
%     G(i,1) = cf(1)*exp(-(((x2_ad(i,1)-cf(2)).^2)/((cf(3)).^2)));
%     L(i,1) = cf(4)*(cf(5)/((x2_ad(i,1)-cf(2)).^2+((cf(5)).^2)));
%     FD(i,1) = ((1+ exp(x2_ad(i,1).*1000/(25+cf(6)))).^(-1));
%     BG(i,1) = cf(7)*Tgd_ads(i,1);
% %     vgt(i,1) = conv(G(i,1),L(i,1),"same")*FD(i,1);
% 
% end
% fitfn = conv(G,L,"same").*FD+BG;
% % fitfn2 = vgt;
% % q = cumtrapz(vgt);
% FW_v = 0.5346*cf(5) +sqrt(0.2166*(cf(5)).^2 +(cf(3)).^2); %FWHM of voigt profile from literature
% HW_v = 0.5*FW_v; %HWHM of voigt profile
% C_v = cf(2); %center of Voigt profile from literature
% [BE_v,loc]= max(fitfn(:,1));
% Mx_v = x2_ad(loc,1); %maxima position of voigt function
% param_HCN = [C_v HW_v T_eff]; %save Voigt parameters to an array
% 
% figure
% plot(x2_ad,y2_ad,'ob',x2_ad,fitfn,'-r','LineWidth',2)
% legend('0.4 ML HATCN/22ML Ag','fit','Location','northwest')
% xlabel('Binding Energy(eV)')
% ylabel('Counts(a.u.)')
% xlim([-0.3 0.3])
% ylim([0 1.2])
% % xline(pk1,'-r')
% a = get(gca,'XTickLabel'); 
% set(gca,'fontsize',18,'FontWeight','bold');
% b = get(gca,'XTick');
% set(gca,'LineWidth',2);
% c = get(gca,'legend');
% set(gca,'fontsize',14);
% legend('boxoff')
% box on