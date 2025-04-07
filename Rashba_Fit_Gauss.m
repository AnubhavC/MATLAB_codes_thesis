% Fitting the Rashba SS

data = load('180HzBiAg2_intp10.mat'); %substrate data
% data = load('HATCN_BiAg2_1114_intp2.mat');
data = data.cube;
k = data(:,:,2);
BE = data(:,:,3);
I = data(:,:,4);

cols = length(k(1,:));
rows = length(k(:,1));
kpar = k(end,:);
BE_range = BE(:,1);

figure  
surf(k,BE,I)
caxis([1400 3500])
colormap('jet')
shading flat
brighten(0);
xlim([k(end,1) k(end,end)])
ylim([BE(1,1) BE(end,1)])
view(0,89.9)
set(gcf,'position',[100,100,550,600])
xlabel('K_|_|(Å^-^1)')
ylabel('BE(eV)')
xline(0,'g')
yline(0,'g')
hold off

%select range of BE and Kpar for band to be fit
BE_start = 1;
BE_mid = 124;
BE_end = length(BE(:,1));
K_start = 1;
K_gamma = find(k(1,:)==0);
K_end = 121;
BE_range = BE(:,1);
BE_low = BE(BE_mid:BE_end,1);
BE_high = BE(BE_start:BE_mid,1);

%Tougaard BG
%~~~~Tougaard BG for substrate
Tgd_sub = zeros(length(BE(:,1)),length(BE(1,:)));
I_sub = zeros(length(BE(:,1)),length(BE(1,:)));
dE = abs(BE_range(1,1)-BE_range(2,1));
P = 20; %Tougaard parameter
for j=1:K_end
    for p = 1:length(BE_range(:,1))
        for q = p:length(BE_range(:,1))
            Tgd_sub(p,j) = Tgd_sub(p,j) + I(q,j)*dE*((BE_range(q,1)-BE_range(p,1))/(P+(BE_range(q,1)-BE_range(p,1).^2).^2));
        end
    end
    B1 = I(1,j)/Tgd_sub(1,j);
    Tgd_sub(:,j) = B1*Tgd_sub(:,j);
    I_sub(:,j) = I(:,j)-Tgd_sub(:,j);% BG subtracted intensities
end

I_smooth = zeros(length(BE(:,1)),length(BE(1,:)));
for j=1:K_end
    I_smooth(:,j) = smoothdata(I(:,j),"sgolay",1);

end

% figure  
% surf(k,BE,I_smooth)
% caxis([700 1400])
% colormap('hot')
% shading flat
% brighten(0);
% xlim([k(end,1) k(end,end)])
% ylim([BE(1,1) BE(end,1)])
% view(0,89.9)
% set(gcf,'position',[100,100,550,600])
% xlabel('K_|_|(Å^-^1)')
% ylabel('BE(eV)')
% xline(0,'g')
% yline(0,'g')
% hold off

%~~~~Find the peak positions at each k||~~~~~

%~~~~~method 1~~~~~~
%divide the BE space into 2 regions: from EF to around -0.55 eV, and ~ -0.55 eV to
%-2 eV..fit each peak separately

%Tougaard BG for low BE region
I_low = I_smooth(BE_mid:BE_end,K_gamma:K_end);
Tgd_low = zeros(length(BE(BE_mid:BE_end,1)),length(BE(1,K_gamma:K_end)));

dE = abs(BE_range(1,1)-BE_range(2,1));
P = 20; %Tougaard parameter
for j=1:K_end-K_gamma
    for p = 1:length(BE_range(BE_mid:BE_end,1))
        for q = p:length(BE_range(BE_mid:BE_end,1))
            Tgd_low(p,j) = Tgd_low(p,j) + I_low(q,j)*dE*((BE_range(q+BE_mid-1,1)-BE_range(p+BE_mid-1,1))/(P+(BE_range(q+BE_mid-1,1)-BE_range(p+BE_mid-1,1).^2).^2));
        end
    end
    B_low = I_low(1,j)/Tgd_low(1,j);
    Tgd_low(:,j) = B_low*Tgd_low(:,j);
    
end
% plot(BE_range(BE_mid:BE_end,1),I_smooth(BE_mid:BE_end,18),BE_range(BE_mid:BE_end,1),Tgd_low(:,5))
% for r=1:length(Tgd_low(1,:))
%     plot(BE_range(BE_mid:BE_end,1),I_smooth(BE_mid:BE_end,K_gamma+r-1),BE_range(BE_mid:BE_end,1),Tgd_low(:,r))
%     hold on
% end
% hold off
%fit parameters
pk = zeros(1,length(I_low(1,:)));%peak BE
HW = zeros(1,length(I_low(1,:)));%HWHM
R2 = zeros(1,length(I_low(1,:)));%R^2
fitfn = zeros(length(BE_low(:,1)),length(I_low(1,:))); %fitted output
FD = zeros(length(BE_low(:,1)),length(I_low(1,:))); % F-D distribution
SS = zeros(length(BE_low(:,1)),length(I_low(1,:))); % fitted Gaussian for Rashba SS

for q=1:length(I_low(1,:))
    BE_fit = BE(BE_mid:BE_end,K_start);

    I_fit = I_low(:,q);
    [mx,loc] = max(I_fit(:,1)); %gaussian peak should be close to the intensity maxima
    ctr = BE_fit(loc,1); %BE for intensity maxima at a given kpar
    tol = 0.05; %tolerance(in eV) for deviation from intensity maxima
    Tgd_fit = Tgd_low(:,q);


    start_param = [1000,ctr,0.1,0.5,I_fit(1,1)];
    lower_param = [0,ctr-tol,0,0,0];
    upper_param = [max(max(I_fit(:,:))),ctr+tol,0.5,1,I_fit(1,1)];

    ftopt = fitoptions('Method','NonlinearLeastSquares','StartPoint',start_param,'Lower',lower_param,'Upper',upper_param);
    ftmethod = fittype("A*exp(-(BE_fit-b).^2/(c).^2)+ (E*Tgd_fit)+ K*((1 + exp(71.*BE_fit)).^(-1)) ",dependent="I_fit",independent="BE_fit",problem="Tgd_fit",coefficients=["A" "b" "c" "E" "K"],options= ftopt);

    % ftmethod = fittype("(A)*((G.^2)/((BE_fit-b).^2+(G).^2)).*((1 + exp(71.*BE_fit)).^(-1))  + E.*Tgd_fit + (d*BE_fit + K)",dependent="I_fit",independent="BE_fit",problem="Tgd_fit",coefficients=["A" "G" "b" "E" "d" "K"],options= ftopt);
    [f,gof] = fit(BE_fit,I_fit,ftmethod,'problem',Tgd_fit);

    cf = coeffvalues(f);
    pk(1,q) = cf(2);
    HW(1,q) = cf(3);
    R2(1,q) = gof.rsquare;
    
    for i=1:length(BE_fit(:,1))
        fitfn(i,q) = cf(1)*exp(-(BE_fit(i,1)-cf(2)).^2/(cf(3)).^2)+ cf(4)*Tgd_fit(i,1)+ cf(5)*((1 + exp(71.*BE_fit(i,1))).^(-1)) ;
        FD(i,q) = cf(5)*((1 + exp(71.*BE_fit(i,1))).^(-1)) ;
        SS(i,q) = cf(1)*exp(-(BE_fit(i,1)-cf(2)).^2/(cf(3)).^2);
    end
%     plot(BE_fit,SS,'*',BE_fit,I_fit)
%     hold on
       
end

% figure  
% surf(k(BE_mid:BE_end,K_gamma:K_end),BE(BE_mid:BE_end,K_gamma:K_end),fitfn)
% caxis([700 1400])
% colormap('jet')
% shading flat
% brighten(0);
% xlim([k(end,1) k(end,end)])
% ylim([BE(1,1) BE(end,1)])
% view(0,89.9)
% set(gcf,'position',[100,100,550,600])
% xlabel('K_|_|(Å^-^1)')
% ylabel('BE(eV)')
% xline(0,'g')
% yline(0,'g')
% hold off

% figure  
% surf(k(BE_mid:BE_end,K_gamma:K_end),BE(BE_mid:BE_end,K_gamma:K_end),I_smooth(BE_mid:BE_end,K_gamma:K_end))
% caxis([700 1400])
% colormap('parula')
% shading flat
% brighten(0);
% xlim([k(end,1) k(end,end)])
% ylim([BE(1,1) BE(end,1)])
% view(0,89.9)
% set(gcf,'position',[100,100,550,600])
% xlabel('K_|_|(Å^-^1)')
% ylabel('BE(eV)')
% xline(0,'g')
% yline(0,'g')
% hold off

% figure  
% surf(k,BE,I_smooth)
% caxis([700 1100])
% colormap('hot')
% shading flat
% brighten(0);
% xlim([k(end,1) k(end,end)])
% ylim([BE(1,1) BE(end,1)])
% view(0,89.9)
% set(gcf,'position',[100,100,550,600])
% xlabel('K_|_|(Å^-^1)')
% ylabel('BE(eV)')
% xline(0,'g')
% yline(0,'g')

% hold on
% h = surf(k,BE,I_smooth);
% z_max = max(max(get(h,'Zdata')));
% z_min = min(min(get(h,'Zdata')));
% on = ones(length(kpar(:,1)),length(kpar(1,K_gamma:K_end)));
% pk_prime = flip(pk);
% plot3(kpar(1,K_gamma:K_end),pk,z_max*on,'ob',kpar(1,K_start:K_gamma),pk_prime,z_max*on,'ob')
% xline(0.142,'g')
% hold off
    
%~~~~Fitting high BE region~~~~~~~

%Tougaard BG for high BE region
I_high = I_smooth(BE_start:BE_mid,K_gamma:K_end);
Tgd_high = zeros(length(BE(BE_start:BE_mid,1)),length(BE(1,K_gamma:K_end)));

dE = abs(BE_range(1,1)-BE_range(2,1));
P = 20; %Tougaard parameter
for j=1:K_end-K_gamma
    for p = 1:length(BE_range(BE_start:BE_mid,1))
        for q = p:length(BE_range(BE_start:BE_mid,1))
            Tgd_high(p,j) = Tgd_high(p,j) + I_high(q,j)*dE*((BE_range(q+BE_start-1,1)-BE_range(p+BE_start-1,1))/(P+(BE_range(q+BE_start-1,1)-BE_range(p+BE_start-1,1).^2).^2));
        end
    end
    B_high = I_high(1,j)/Tgd_high(1,j);
    Tgd_high(:,j) = B_high*Tgd_high(:,j);
    
end
% plot(BE_range(BE_mid:BE_end,1),I_smooth(BE_mid:BE_end,18),BE_range(BE_mid:BE_end,1),Tgd_low(:,5))
% for r=1:length(Tgd_high(1,:))
%     plot(BE_range(BE_start:BE_mid,1),I_smooth(BE_start:BE_mid,K_gamma+r-1),BE_range(BE_start:BE_mid,1),Tgd_high(:,r))
%     hold on
% end
hold off

%fit parameters
pk_high = zeros(1,length(I_high(1,:)));%peak BE
HW_high = zeros(1,length(I_high(1,:)));%HWHM
R2_high = zeros(1,length(I_high(1,:)));%R^2
fitfn_high = zeros(length(BE_high(:,1)),length(I_high(1,:))); %fitted output
FD_high = zeros(length(BE_high(:,1)),length(I_high(1,:))); % F-D distribution
SS_high = zeros(length(BE_high(:,1)),length(I_high(1,:))); % fitted Gaussian for Rashba SS

for q=1:length(I_high(1,:))
    BE_fit_high = BE(BE_start:BE_mid,K_start);

    I_fit_high = I_high(:,q);
    [mx,loc] = max(I_fit_high(:,1)); %gaussian peak should be close to the intensity maxima
    ctr = BE_fit_high(loc,1); %BE for intensity maxima at a given kpar
    tol = 0.05; %tolerance(in eV) for deviation from intensity maxima 
    Tgd_fit_high = Tgd_high(:,q);
    counter = zeros(1,length(I_high(1,:)));

%     if max(Tgd_fit_high(:,1))>= mx
%         Tgd_fit_high = Tgd_high(:,1);
%         counter(1,q) =1;
%     end

    start_param = [1000,ctr,0.1,0.5,I_fit(1,1)];
    lower_param = [0,ctr-tol,0,0,0];
    upper_param = [max(max(I_fit(:,:))),ctr+tol,0.5,1,I_fit(1,1)];

    ftopt = fitoptions('Method','NonlinearLeastSquares','StartPoint',start_param,'Lower',lower_param,'Upper',upper_param);
    ftmethod = fittype("A*exp(-(BE_fit_high-b).^2/(c).^2)+ (E*Tgd_fit_high)+ K",dependent="I_fit_high",independent="BE_fit_high",problem="Tgd_fit_high",coefficients=["A" "b" "c" "E" "K"],options= ftopt);

    % ftmethod = fittype("(A)*((G.^2)/((BE_fit-b).^2+(G).^2)).*((1 + exp(71.*BE_fit)).^(-1))  + E.*Tgd_fit + (d*BE_fit + K)",dependent="I_fit",independent="BE_fit",problem="Tgd_fit",coefficients=["A" "G" "b" "E" "d" "K"],options= ftopt);
    [f,gof] = fit(BE_fit_high,I_fit_high,ftmethod,'problem',Tgd_fit_high);

    cf = coeffvalues(f);
    pk_high(1,q) = cf(2);
    HW_high(1,q) = cf(3);
    R2_high(1,q) = gof.rsquare;
    
    for i=1:length(BE_fit_high(:,1))
        fitfn_high(i,q) = cf(1)*exp(-(BE_fit_high(i,1)-cf(2)).^2/(cf(3)).^2)+ cf(4)*Tgd_fit_high(i,1)+ cf(5);
        SS_high(i,q) = cf(1)*exp(-(BE_fit_high(i,1)-cf(2)).^2/(cf(3)).^2);
    end
    
   
end
%     l=length(I_high(1,:));
%     for w=l-2:l
%         plot(BE_fit_high,fitfn_high(:,w))
%         hold on
%     end
pk_prime_high = flip(pk_high); %array for negative kpar
pk_gamma = 0.5*(pk_high(1,1)+pk(1,1)); %peak position for Gamma
prbl = zeros(1,length(kpar(1,K_start:K_end)));

for q=1:length(kpar(1,K_start:K_end))
    if q<K_gamma
        prbl(1,q) = pk_prime_high(1,q);
        
    elseif q>K_gamma  
        prbl(1,q) = pk(1,q-K_gamma);
        
    else
        prbl(1,q) = pk_gamma;
        
    end
end
fitfn_prime = fliplr(fitfn);
fitfn_high_prime = fliplr(fitfn_high);


% plot(BE_range,I_smooth(:,K_gamma))

% figure  
% surf(k,BE,I_smooth)
% caxis([700 1100])
% colormap('hot')
% shading flat
% brighten(0);
% xlim([k(end,1) k(end,end)])
% ylim([BE(1,1) BE(end,1)])
% view(0,89.9)
% set(gcf,'position',[100,100,550,600])
% xlabel('K_|_|(Å^-^1)')
% ylabel('BE(eV)')
% xline(0,'g')
% yline(0,'g')
% 
% hold on
% h = surf(k,BE,I_smooth);
% z_max = max(max(get(h,'Zdata')));
% z_min = min(min(get(h,'Zdata')));
% on = ones(length(kpar(:,1)),length(kpar(1,K_start:K_end)));
% 
% plot3(kpar(1,K_start:K_end),prbl,z_max*on,'ob')
% hold off

%Fit using parabola to find Rashba parameters

%parabolae
h_bar =9;
m_eff = -0.33;
k_offset = 0.12;
BE_offset = 0.1;

k_fit = transpose(kpar(1,K_start:K_end));
prbl = transpose(prbl);

% y1 = zeros(length(k_fit(:,1)),1);
% y1 = ((h_bar.^2)*((k_fit+k_offset).^2)/2*m_eff)-BE_offset;
% % plot(k_fit,y1,'ob')
% 
% y2 = zeros(length(k_fit(:,1)),1);
% y2 = ((h_bar.^2)*((k_fit-k_offset).^2)/2*m_eff)-BE_offset;
% plot(k_fit,y1,'-ob',k_fit,y2,'-or')
% hold off
% scatter(k_fit,prbl)

%intial guess values
start_param = [((h_bar.^2)/2*m_eff),k_offset,BE_offset];
lower_param = [-Inf,0,-1];
upper_param = [0,k_fit(end,1),1];
tf = excludedata(k_fit,prbl,'indices',K_gamma);
% scatter(k_fit(~tf),prbl(~tf))
% hold on
% plot(k_fit,prbl,'*')

ftopt = fitoptions('Method','NonlinearLeastSquares','Lower',lower_param,'Upper',upper_param,'StartPoint',start_param,'TolX',1e-6,'TolFun',1e-6,'MaxIter',400);
% ftmethod1 = fittype('(M*((x+k).^2))-E','independent',{'x'},'dependent',{'prbol1'},'coefficients',{'M','k','E'},'options',ftopt);
% [f1,gof1] = fit(k_fit(1:K_end,1),prbol1(1:K_end,1),ftmethod1); %fit BE vs k w/ parabola

ftmethod2 = fittype('(M*((k_fit-C).^2))-E','independent',{'k_fit'},'dependent',{'prbl'},'coefficients',{'M','C','E'},'options',ftopt);
[f2,gof2] = fit(k_fit,prbl,ftmethod2,'Exclude',tf); %fit BE vs k w/ parabola

r2 = gof2.rsquare;
cf2 = coeffvalues(f2);
fitfn2 =zeros();
fitfn3 =zeros();

for i=1:length(k_fit(:,1))
    
    fitfn2(i,1) = (cf2(1)*((k_fit(i,1)-cf2(2)).^2))-cf2(3);
    fitfn3(i,1) = (cf2(1)*((k_fit(i,1)+cf2(2)).^2))-cf2(3);

end
%fitted Rashba parameters
h_bar = 6.58*1E-16;
me = 9.1*1E-31;
Rashba_meff = 3.81/cf2(1); %effective mass of band
Rashba_K0 = cf2(2); %parabola center
Rashba_BE = cf2(3); %BE offset
Rashba_crossing = fitfn2(K_gamma,1);

% E_Rashba = abs(((cf2(1)*((Rashba_K0-cf2(2)).^2))-cf2(3))-((cf2(1)*((k_fit(K_gamma,1)-cf2(2)).^2))-cf2(3)));
E_Rashba = ((cf2(1)*((Rashba_K0-cf2(2)).^2))+cf2(3))-((cf2(1)*((Rashba_K0+cf2(2)).^2))+cf2(3));

alphaR = E_Rashba/(2*Rashba_K0);

% plot(k_fit,prbl,'ob',k_fit,fitfn2,'MarkerSize',5,'LineWidth',1)
close all

figure  
surf(k,BE,I_smooth)
caxis([1400 3500])
colormap('hot')
shading flat
brighten(0);
xlim([k(end,1) k(end,K_end)])
ylim([BE(1,1) 0.2])
view(0,89.9)
set(gcf,'position',[100,100,500,605])
xlabel('k_|_|(Å^-^1)','FontSize',16,'Units','Normalized','Position',[0.6 -0.06],'FontWeight','bold')
a = get(gca,'XTickLabel');
set(gca,'fontsize',14,'FontWeight','bold');
b = get(gca,'XTick');
set(gca,'LineWidth',2);
ylabel('Binding Energy(eV)','FontSize',16,'Units','Normalized')
xline(0,'b')
yline(0,'r')
% xline(0.142,'--g')
% xline(Rashba_K0,'g')
% yline(Rashba_crossing,'g')
hold on
% for i=1:length(k_fit(:,1))
%     
%     x_i = k_fit(i,1);
%     xline(x_i,'-g')
%     hold on
% end
h = surf(k,BE,I_smooth);
z_max = max(max(get(h,'Zdata')));
z_min = min(min(get(h,'Zdata')));
on = ones(length(k_fit(:,1)),length(k_fit(1,:)));

plot3(k_fit,fitfn2,z_max*on,'.w',k_fit,fitfn3,z_max*on,'.w')
xlim([k(end,1) k(end,K_end)])
shading flat
hold off

% figure %Fitted ARPES with parabolae  
% surf(k(BE_start:BE_mid,K_gamma:K_end),BE(BE_start:BE_mid,K_gamma:K_end),fitfn_high)
% hold on
% surf(k(BE_mid:BE_end,K_gamma:K_end),BE(BE_mid:BE_end,K_gamma:K_end),fitfn)
% hold on
% surf(k(BE_start:BE_mid,K_start:K_gamma),BE(BE_start:BE_mid,K_start:K_gamma),fitfn_high_prime)
% hold on
% surf(k(BE_mid:BE_end,K_start:K_gamma),BE(BE_mid:BE_end,K_start:K_gamma),fitfn_prime)
% 
% 
% caxis([700 1100])
% colormap('hot')
% shading flat
% brighten(0);
% xlim([k(end,1) k(end,end)])
% ylim([BE(1,1) BE(end,1)])
% view(0,89.9)
% set(gcf,'position',[100,100,550,600])
% xlabel('K_|_|(Å^-^1)')
% ylabel('BE(eV)')
% xline(0,'g')
% yline(0,'g')
% hold on
% 
% h = surf(k(BE_mid:BE_end,K_gamma:K_end),BE(BE_mid:BE_end,K_gamma:K_end),fitfn);
% z_max = max(max(get(h,'Zdata')));
% z_min = min(min(get(h,'Zdata')));
% on = ones(length(k_fit(:,1)),length(k_fit(1,:)));
% 
% plot3(k_fit,fitfn2,z_max*on,'ob',k_fit,fitfn3,z_max*on,'og')
% 
% shading flat
% hold off

%Save fitted parabolae data and RSS parameters
fitfn_BiAg_R = fitfn2;
fitfn_BiAg_L = fitfn3;
RSS_BiAg2 = zeros();
RSS_BiAg2(1,1) =Rashba_meff;
RSS_BiAg2(1,2) =Rashba_K0;
RSS_BiAg2(1,3) =Rashba_BE;%BE offset
RSS_BiAg2(1,4) =abs(Rashba_crossing)-Rashba_BE;%ER
RSS_BiAg2(1,5) = alphaR;

%&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&

data = load('NPT_BiAg2_21_07_24_intp2.mat'); %adsorbate data
% data = load('HATCN_BiAg2_1114_intp2.mat');
data = data.cube;
k = data(:,:,2);
BE = data(:,:,3);
I = data(:,:,4);

cols = length(k(1,:));
rows = length(k(:,1));
kpar = k(end,:);
BE_range = BE(:,1);

% figure  
% surf(k,BE,I)
% caxis([1400 2700])
% colormap('gray')
% shading flat
% brighten(0);
% xlim([k(end,1) k(end,end)])
% ylim([BE(1,1) BE(end,1)])
% view(0,89.9)
% set(gcf,'position',[100,100,550,600])
% xlabel('K_|_|(Å^-^1)')
% ylabel('BE(eV)')
% xline(0,'g')
% yline(0,'g')
% hold on

%select range of BE and Kpar for band to be fit
BE_start = 1;
BE_mid = 115;
BE_end = length(BE(:,1));
K_start = 1;
K_gamma = find(k(1,:)==0);
K_end = 27;
BE_range = BE(:,1);
BE_low = BE(BE_mid:BE_end,1);
BE_high = BE(BE_start:BE_mid,1);

%Tougaard BG
%~~~~Tougaard BG for substrate
Tgd_sub = zeros(length(BE(:,1)),length(BE(1,:)));
I_sub = zeros(length(BE(:,1)),length(BE(1,:)));
dE = abs(BE_range(1,1)-BE_range(2,1));
P = 1643; %Tougaard parameter
for j=1:K_end
    for p = 1:length(BE_range(:,1))
        for q = p:length(BE_range(:,1))
            Tgd_sub(p,j) = Tgd_sub(p,j) + I(q,j)*dE*((BE_range(q,1)-BE_range(p,1))/(P+(BE_range(q,1)-BE_range(p,1).^2).^2));
        end
    end
    B1 = I(1,j)/Tgd_sub(1,j);
    Tgd_sub(:,j) = B1*Tgd_sub(:,j);
    I_sub(:,j) = I(:,j)-Tgd_sub(:,j);% BG subtracted intensities
end

I_smooth = zeros(length(BE(:,1)),length(BE(1,:)));
for j=1:K_end
    I_smooth(:,j) = smoothdata(I(:,j),"sgolay",10);

end

figure  
surf(k,BE,I_smooth)
caxis([1300 2700])
colormap('gray')
shading flat
brighten(0);
xlim([k(end,1) k(end,end)])
ylim([BE(1,1) BE(end,1)])
view(0,89.9)
set(gcf,'position',[100,100,550,600])
xlabel('K_|_|(Å^-^1)')
ylabel('BE(eV)')
xline(0,'g')
yline(0,'g')
hold on
plot3(k_fit,fitfn2,z_max*on,'ob',k_fit,fitfn3,z_max*on,'og')
shading flat
hold off

%~~~~Find the peak positions at each k||~~~~~

%~~~~~method 1~~~~~~
%divide the BE space into 2 regions: from EF to around -0.55 eV, and ~ -0.55 eV to
%-2 eV..fit each peak separately

%Tougaard BG for low BE region
I_low = I_smooth(BE_mid:BE_end,K_gamma:K_end);
Tgd_low = zeros(length(BE(BE_mid:BE_end,1)),length(BE(1,K_gamma:K_end)));

dE = abs(BE_range(1,1)-BE_range(2,1));
P = 1643; %Tougaard parameter
for j=1:K_end-K_gamma
    for p = 1:length(BE_range(BE_mid:BE_end,1))
        for q = p:length(BE_range(BE_mid:BE_end,1))
            Tgd_low(p,j) = Tgd_low(p,j) + I_low(q,j)*dE*((BE_range(q+BE_mid-1,1)-BE_range(p+BE_mid-1,1))/(P+(BE_range(q+BE_mid-1,1)-BE_range(p+BE_mid-1,1).^2).^2));
        end
    end
    B_low = I_low(1,j)/Tgd_low(1,j);
    Tgd_low(:,j) = B_low*Tgd_low(:,j);
    
end
% plot(BE_range(BE_mid:BE_end,1),I_smooth(BE_mid:BE_end,18),BE_range(BE_mid:BE_end,1),Tgd_low(:,5))
% for r=1:length(Tgd_low(1,:))
%     plot(BE_range(BE_mid:BE_end,1),I_smooth(BE_mid:BE_end,K_gamma+r-1),BE_range(BE_mid:BE_end,1),Tgd_low(:,r))
%     hold on
% end
% hold off
%fit parameters
pk = zeros(1,length(I_low(1,:)));%peak BE
HW = zeros(1,length(I_low(1,:)));%HWHM
R2 = zeros(1,length(I_low(1,:)));%R^2
fitfn = zeros(length(BE_low(:,1)),length(I_low(1,:))); %fitted output
FD = zeros(length(BE_low(:,1)),length(I_low(1,:))); % F-D distribution
SS = zeros(length(BE_low(:,1)),length(I_low(1,:))); % fitted Gaussian for Rashba SS

for q=1:length(I_low(1,:))
    BE_fit = BE(BE_mid:BE_end,K_start);

    I_fit = I_low(:,q);
    [mx,loc] = max(I_fit(:,1)); %gaussian peak should be close to the intensity maxima
    ctr = BE_fit(loc,1); %BE for intensity maxima at a given kpar
    tol = 0.05; %tolerance(in eV) for deviation from intensity maxima
    Tgd_fit = Tgd_low(:,q);


    start_param = [1000,ctr,0.1,0.5,I_fit(1,1)];
    lower_param = [0,ctr-tol,0,0,0];
    upper_param = [max(max(I_fit(:,:))),ctr+tol,0.5,1,I_fit(1,1)];

    ftopt = fitoptions('Method','NonlinearLeastSquares','StartPoint',start_param,'Lower',lower_param,'Upper',upper_param);
    ftmethod = fittype("A*exp(-(BE_fit-b).^2/(c).^2)+ (E*Tgd_fit)+ K*((1 + exp(71.*BE_fit)).^(-1)) ",dependent="I_fit",independent="BE_fit",problem="Tgd_fit",coefficients=["A" "b" "c" "E" "K"],options= ftopt);

    % ftmethod = fittype("(A)*((G.^2)/((BE_fit-b).^2+(G).^2)).*((1 + exp(71.*BE_fit)).^(-1))  + E.*Tgd_fit + (d*BE_fit + K)",dependent="I_fit",independent="BE_fit",problem="Tgd_fit",coefficients=["A" "G" "b" "E" "d" "K"],options= ftopt);
    [f,gof] = fit(BE_fit,I_fit,ftmethod,'problem',Tgd_fit);

    cf = coeffvalues(f);
    pk(1,q) = cf(2);
    HW(1,q) = cf(3);
    R2(1,q) = gof.rsquare;
    
    for i=1:length(BE_fit(:,1))
        fitfn(i,q) = cf(1)*exp(-(BE_fit(i,1)-cf(2)).^2/(cf(3)).^2)+ cf(4)*Tgd_fit(i,1)+ cf(5)*((1 + exp(71.*BE_fit(i,1))).^(-1)) ;
        FD(i,q) = cf(5)*((1 + exp(71.*BE_fit(i,1))).^(-1)) ;
        SS(i,q) = cf(1)*exp(-(BE_fit(i,1)-cf(2)).^2/(cf(3)).^2);
    end
%     plot(BE_fit,SS,'*',BE_fit,I_fit)
%     hold on
       
end

% figure  
% surf(k(BE_mid:BE_end,K_gamma:K_end),BE(BE_mid:BE_end,K_gamma:K_end),fitfn)
% caxis([700 1400])
% colormap('jet')
% shading flat
% brighten(0);
% xlim([k(end,1) k(end,end)])
% ylim([BE(1,1) BE(end,1)])
% view(0,89.9)
% set(gcf,'position',[100,100,550,600])
% xlabel('K_|_|(Å^-^1)')
% ylabel('BE(eV)')
% xline(0,'g')
% yline(0,'g')
% hold off

% figure  
% surf(k(BE_mid:BE_end,K_gamma:K_end),BE(BE_mid:BE_end,K_gamma:K_end),I_smooth(BE_mid:BE_end,K_gamma:K_end))
% caxis([700 1400])
% colormap('parula')
% shading flat
% brighten(0);
% xlim([k(end,1) k(end,end)])
% ylim([BE(1,1) BE(end,1)])
% view(0,89.9)
% set(gcf,'position',[100,100,550,600])
% xlabel('K_|_|(Å^-^1)')
% ylabel('BE(eV)')
% xline(0,'g')
% yline(0,'g')
% hold off

% figure  
% surf(k,BE,I_smooth)
% caxis([700 1100])
% colormap('hot')
% shading flat
% brighten(0);
% xlim([k(end,1) k(end,end)])
% ylim([BE(1,1) BE(end,1)])
% view(0,89.9)
% set(gcf,'position',[100,100,550,600])
% xlabel('K_|_|(Å^-^1)')
% ylabel('BE(eV)')
% xline(0,'g')
% yline(0,'g')

% hold on
% h = surf(k,BE,I_smooth);
% z_max = max(max(get(h,'Zdata')));
% z_min = min(min(get(h,'Zdata')));
% on = ones(length(kpar(:,1)),length(kpar(1,K_gamma:K_end)));
% pk_prime = flip(pk);
% plot3(kpar(1,K_gamma:K_end),pk,z_max*on,'ob',kpar(1,K_start:K_gamma),pk_prime,z_max*on,'ob')
% xline(0.142,'g')
% hold off
    
%~~~~Fitting high BE region~~~~~~~

%Tougaard BG for high BE region
I_high = I_smooth(BE_start:BE_mid,K_gamma:K_end);
Tgd_high = zeros(length(BE(BE_start:BE_mid,1)),length(BE(1,K_gamma:K_end)));

dE = abs(BE_range(1,1)-BE_range(2,1));
P = 1643; %Tougaard parameter
for j=1:K_end-K_gamma
    for p = 1:length(BE_range(BE_start:BE_mid,1))
        for q = p:length(BE_range(BE_start:BE_mid,1))
            Tgd_high(p,j) = Tgd_high(p,j) + I_high(q,j)*dE*((BE_range(q+BE_start-1,1)-BE_range(p+BE_start-1,1))/(P+(BE_range(q+BE_start-1,1)-BE_range(p+BE_start-1,1).^2).^2));
        end
    end
    B_high = I_high(1,j)/Tgd_high(1,j);
    Tgd_high(:,j) = B_high*Tgd_high(:,j);
    
end
% plot(BE_range(BE_mid:BE_end,1),I_smooth(BE_mid:BE_end,18),BE_range(BE_mid:BE_end,1),Tgd_low(:,5))
% for r=1:length(Tgd_high(1,:))
%     plot(BE_range(BE_start:BE_mid,1),I_smooth(BE_start:BE_mid,K_gamma+r-1),BE_range(BE_start:BE_mid,1),Tgd_high(:,r))
%     hold on
% end
hold off

%fit parameters
pk_high = zeros(1,length(I_high(1,:)));%peak BE
HW_high = zeros(1,length(I_high(1,:)));%HWHM
R2_high = zeros(1,length(I_high(1,:)));%R^2
fitfn_high = zeros(length(BE_high(:,1)),length(I_high(1,:))); %fitted output
FD_high = zeros(length(BE_high(:,1)),length(I_high(1,:))); % F-D distribution
SS_high = zeros(length(BE_high(:,1)),length(I_high(1,:))); % fitted Gaussian for Rashba SS

for q=1:length(I_high(1,:))
    BE_fit_high = BE(BE_start:BE_mid,K_start);

    I_fit_high = I_high(:,q);
    [mx,loc] = max(I_fit_high(:,1)); %gaussian peak should be close to the intensity maxima
    ctr = BE_fit_high(loc,1); %BE for intensity maxima at a given kpar
    tol = 0.05; %tolerance(in eV) for deviation from intensity maxima 
    Tgd_fit_high = Tgd_high(:,q);
    counter = zeros(1,length(I_high(1,:)));

%     if max(Tgd_fit_high(:,1))>= mx
%         Tgd_fit_high = Tgd_high(:,1);
%         counter(1,q) =1;
%     end

    start_param = [1000,ctr,0.1,0.5,I_fit(1,1)];
    lower_param = [0,ctr-tol,0,0,0];
    upper_param = [max(max(I_fit(:,:))),ctr+tol,0.5,1,I_fit(1,1)];

    ftopt = fitoptions('Method','NonlinearLeastSquares','StartPoint',start_param,'Lower',lower_param,'Upper',upper_param);
    ftmethod = fittype("A*exp(-(BE_fit_high-b).^2/(c).^2)+ (E*Tgd_fit_high)+ K",dependent="I_fit_high",independent="BE_fit_high",problem="Tgd_fit_high",coefficients=["A" "b" "c" "E" "K"],options= ftopt);

    % ftmethod = fittype("(A)*((G.^2)/((BE_fit-b).^2+(G).^2)).*((1 + exp(71.*BE_fit)).^(-1))  + E.*Tgd_fit + (d*BE_fit + K)",dependent="I_fit",independent="BE_fit",problem="Tgd_fit",coefficients=["A" "G" "b" "E" "d" "K"],options= ftopt);
    [f,gof] = fit(BE_fit_high,I_fit_high,ftmethod,'problem',Tgd_fit_high);

    cf = coeffvalues(f);
    pk_high(1,q) = cf(2);
    HW_high(1,q) = cf(3);
    R2_high(1,q) = gof.rsquare;
    
    for i=1:length(BE_fit_high(:,1))
        fitfn_high(i,q) = cf(1)*exp(-(BE_fit_high(i,1)-cf(2)).^2/(cf(3)).^2)+ cf(4)*Tgd_fit_high(i,1)+ cf(5);
        SS_high(i,q) = cf(1)*exp(-(BE_fit_high(i,1)-cf(2)).^2/(cf(3)).^2);
    end
    
   
end
%     l=length(I_high(1,:));
%     for w=l-2:l
%         plot(BE_fit_high,fitfn_high(:,w))
%         hold on
%     end
pk_prime_high = flip(pk_high); %array for negative kpar
pk_gamma = 0.5*(pk_high(1,1)+pk(1,1)); %peak position for Gamma
prbl = zeros(1,length(kpar(1,K_start:K_end)));

for q=1:length(kpar(1,K_start:K_end))
    if q<K_gamma
        prbl(1,q) = pk_prime_high(1,q);
        
    elseif q>K_gamma  
        prbl(1,q) = pk(1,q-K_gamma);
        
    else
        prbl(1,q) = pk_gamma;
        
    end
end
fitfn_prime = fliplr(fitfn);
fitfn_high_prime = fliplr(fitfn_high);


% plot(BE_range,I_smooth(:,K_gamma))

% figure  
% surf(k,BE,I_smooth)
% caxis([700 1100])
% colormap('hot')
% shading flat
% brighten(0);
% xlim([k(end,1) k(end,end)])
% ylim([BE(1,1) BE(end,1)])
% view(0,89.9)
% set(gcf,'position',[100,100,550,600])
% xlabel('K_|_|(Å^-^1)')
% ylabel('BE(eV)')
% xline(0,'g')
% yline(0,'g')
% 
% hold on
% h = surf(k,BE,I_smooth);
% z_max = max(max(get(h,'Zdata')));
% z_min = min(min(get(h,'Zdata')));
% on = ones(length(kpar(:,1)),length(kpar(1,K_start:K_end)));
% 
% plot3(kpar(1,K_start:K_end),prbl,z_max*on,'ob')
% hold off

%Fit using parabola to find Rashba parameters

%parabolae
h_bar =9;
m_eff = -0.3;
k_offset = 0.13;
BE_offset = 0.1;

k_fit = transpose(kpar(1,K_start:K_end));
prbl = transpose(prbl);

% y1 = zeros(length(k_fit(:,1)),1);
% y1 = ((h_bar.^2)*((k_fit+k_offset).^2)/2*m_eff)-BE_offset;
% % plot(k_fit,y1,'ob')
% 
% y2 = zeros(length(k_fit(:,1)),1);
% y2 = ((h_bar.^2)*((k_fit-k_offset).^2)/2*m_eff)-BE_offset;
% plot(k_fit,y1,'-ob',k_fit,y2,'-or')
% hold off
% scatter(k_fit,prbl)

%intial guess values
start_param = [((h_bar.^2)/2*m_eff),k_offset,BE_offset];
lower_param = [-Inf,0,-1];
upper_param = [0,k_fit(end,1),1];
tf = excludedata(k_fit,prbl,'indices',K_gamma);
% scatter(k_fit(~tf),prbl(~tf))
% hold on
% plot(k_fit,prbl,'*')

ftopt = fitoptions('Method','NonlinearLeastSquares','Lower',lower_param,'Upper',upper_param,'StartPoint',start_param,'TolX',1e-6,'TolFun',1e-6,'MaxIter',400);
% ftmethod1 = fittype('(M*((x+k).^2))-E','independent',{'x'},'dependent',{'prbol1'},'coefficients',{'M','k','E'},'options',ftopt);
% [f1,gof1] = fit(k_fit(1:K_end,1),prbol1(1:K_end,1),ftmethod1); %fit BE vs k w/ parabola

ftmethod2 = fittype('(M*((k_fit-C).^2))-E','independent',{'k_fit'},'dependent',{'prbl'},'coefficients',{'M','C','E'},'options',ftopt);
[f2,gof2] = fit(k_fit,prbl,ftmethod2,'Exclude',tf); %fit BE vs k w/ parabola

r2 = gof2.rsquare;
cf2 = coeffvalues(f2);
fitfn2 =zeros();
fitfn3 =zeros();

for i=1:length(k_fit(:,1))
    
    fitfn2(i,1) = (cf2(1)*((k_fit(i,1)-cf2(2)).^2))-cf2(3);
    fitfn3(i,1) = (cf2(1)*((k_fit(i,1)+cf2(2)).^2))-cf2(3);

end
%fitted Rashba parameters
h_bar = 6.58*1E-16;
me = 9.1*1E-31;
Rashba_meff = 3.81/cf2(1); %effective mass of band
Rashba_K0 = cf2(2); %parabola center
Rashba_BE = cf2(3); %BE offset
Rashba_crossing = fitfn2(K_gamma,1);

% E_Rashba = abs(((cf2(1)*((Rashba_K0-cf2(2)).^2))-cf2(3))-((cf2(1)*((k_fit(K_gamma,1)-cf2(2)).^2))-cf2(3)));
E_Rashba = ((cf2(1)*((Rashba_K0-cf2(2)).^2))+cf2(3))-((cf2(1)*((Rashba_K0+cf2(2)).^2))+cf2(3));

alphaR = E_Rashba/(2*Rashba_K0);

% plot(k_fit,prbl,'ob',k_fit,fitfn2,'MarkerSize',5,'LineWidth',1)
close all

figure  
surf(k,BE,I_smooth)
caxis([1400 2700])
colormap('hot')
shading flat
brighten(0);
xlim([k(end,1) k(end,end)])
ylim([BE(1,1) BE(end,1)])
view(0,89.9)
set(gcf,'position',[100,100,550,600])
xlabel('K_|_|(Å^-^1)')
ylabel('BE(eV)')
xline(0,'g')
yline(0,'g')
% xline(0.142,'--g')
% xline(Rashba_K0,'g')
% yline(Rashba_crossing,'g')
hold on
for i=1:length(k_fit(:,1))
    
    x_i = k_fit(i,1);
    xline(x_i,'-g')
    hold on
end
h = surf(k,BE,I_smooth);
z_max = max(max(get(h,'Zdata')));
z_min = min(min(get(h,'Zdata')));
on = ones(length(k_fit(:,1)),length(k_fit(1,:)));

plot3(k_fit,fitfn2,z_max*on,'ob',k_fit,fitfn3,z_max*on,'og')
shading flat
hold off

%Save fitted parabolae data
fitfn_HCNBiAg_R = fitfn2;
fitfn_HCNBiAg_L = fitfn3;
RSS_HCNBiAg2 = zeros();
RSS_HCNBiAg2(1,1) =Rashba_meff;
RSS_HCNBiAg2(1,2) =Rashba_K0;
RSS_HCNBiAg2(1,3) =Rashba_BE;%BE offset
RSS_HCNBiAg2(1,4) =abs(Rashba_crossing)-Rashba_BE;%ER
RSS_HCNBiAg2(1,5) = alphaR;

err=0.04*ones(size(fitfn_HCNBiAg_L(:,1)));
% plot(k_fit,fitfn_BiAg_L,'ob',k_fit,fitfn_BiAg_R,'or',k_fit,fitfn_HCNBiAg_L,'*b',k_fit,fitfn_HCNBiAg_R,'*r')
plot(k_fit,fitfn_BiAg_L,'ob')
hold on
plot(k_fit,fitfn_BiAg_R,'or')
hold on
errorbar(k_fit,fitfn_HCNBiAg_L,err,'*b')
hold on
errorbar(k_fit,fitfn_HCNBiAg_R,err,'*r')
hold on

xlabel('K_|_|(Å^-^1)')
ylabel('BE(eV)')
% legend('BiAg_2','BiAg_2','HATCN/BiAg_2','HATCN/BiAg_2','Location','south')
legend('BiAg_2','BiAg_2','NPT/BiAg_2','NPT/BiAg_2','Location','south')


