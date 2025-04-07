%code to import sp band data and normalize data
files = dir('*_SS*');
l = length(files);

BE =zeros();
I =zeros();

BE_start =1;
BE_norm = 1; %choose BE for  normalization

for i=1:l
    A= dlmread(files(i).name);
    x = length(A(BE_start:end,1));
    BE(1:x,i)= A(BE_start:end,1);
    I(1:x,i)= A(BE_start:end,2);
    I(:,i) = I(:,i) - min(I(:,i));
    
   
end
BE = BE-22.08;

I_norm =I;
norm_factor = zeros();
% for i=1:l
%     norm_factor(1,i) = I(BE_norm,3)/I(BE_norm,i);
%     I_norm(:,i) = I(:,i).*norm_factor(1,i);
% end
 
for i=1:l
    
    plot(BE(:,i),I_norm(:,i),'.','MarkerSize',15)
%     plot(BE(:,i),I(:,i))
%     ylim([0,1.25])
    xlim([-2.2,0.5])
    set(gcf,'position',[100,100,550,600])
    xlabel('BE(eV)')
    ylabel('Counts')
    hold on  
end
% Legend = {'1.2 ML','0.6 ML','1 ML','1.9 ML','Bi(110)/Ag','0.3 ML'};
% legend(Legend)
hold off

%fitting NPT peak
x_sub = BE(:,2);
x_NPT = BE(:,4);
y_sub = I_norm(:,2);
y_NPT = I_norm(:,3);

plot(x_sub,y_sub,x_NPT,y_NPT)
xlim([-2,0.4])
set(gcf,'position',[100,100,550,600])
xlabel('BE(eV)')
ylabel('Counts')
legend('Bi110','NPT')

start_param = [1500,-1,0.3,0.5]; 
lower_param = [0,-2,0,0];
upper_param = [max(y_NPT(:,1)),-0.5,1,5];

ftopt = fitoptions('Method','NonlinearLeastSquares','StartPoint',start_param,'Lower',lower_param,'Upper',upper_param);
ftmethod = fittype('A*exp(-(((x_NPT-b).^2)/(c.^2)))+ D*(y_sub)','problem',{'y_sub'},'independent',{'x_NPT'},'dependent',{'y_NPT'},'coefficients',{'A','b','c','D'},'options',ftopt);

[f,gof] = fit(x_NPT,y_NPT,ftmethod,'problem',y_sub);
cf = coeffvalues(f);
r2 = gof.rsquare;
pk = cf(2);
HW = cf(3);
se=confint(f);
err=abs(pk-se(2,2));
for i=1:length(x_NPT(:,1))
    IS(i,1) = cf(1)*exp(-(x_NPT(i,1)-cf(2)).^2/(cf(3)).^2);
    fitfn(i,1) = IS(i,1)+cf(4)*y_sub(i,1);
    
end
plot(x_NPT,y_NPT,'g',x_NPT,fitfn,'-k',x_NPT,IS,'or')
xlim([-2,0.4])
set(gcf,'position',[100,100,550,600])
 set(gcf,'position',[100,100,500,600])
    xlabel('Binding Energy(eV)')
    ylabel('Counts(a.u.)')
    set(gca,'YTick', [])
    a = get(gca,'XTickLabel');
    set(gca,'fontsize',14,'FontWeight','bold');
    b = get(gca,'XTick');
    set(gca,'LineWidth',2,'XMinorTick','on');
    legend('1 ML NO_2-PyT/BiAg_2','fit','ICT',Location='southwest')
    legend("boxoff")