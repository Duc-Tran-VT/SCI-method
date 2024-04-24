% function results = ADV(xin,yin,Cin,xlab,ylab)
clear all; clc; %close all;

load('All.mat'); load('S1S2.mat');

xin = all.adv; x_cal = data.adv; xlab = 'Wetlabs'; 
yin = all.hydro4; y_cal = data.hydro4; ylab = 'Wetlabs'; name = 'ADV_Wetlabs';

Cin = all.C; y_cal_log = 10*log10(y_cal);

fac0 = 0.9*[15,25,50,100,150,200,11.25,18.75,37.5,75,112.5,150,7.5,12.5,25,...
    50,75,100,3.75,6.25,12.5,25,37.5,50,0.001,0.001,0.001,0.001,0.001,0.001,...
        15,25,50,100,150,200,11.25,18.75,37.5,75,112.5,150,7.5,12.5,25,...
    50,75,100,3.75,6.25,12.5,25,37.5,50,0.001,0.001,0.001,0.001,0.001,0.001]'; 

% id = [4 7:11 13:24 28 34 37:45 49:54 58]'; % part 1 from 1:29, part 2 from 30:59 (in new index)
id = [1:11 13:60]'; % part 1 from 1:29, part 2 from 30:59 (in new index)
x=xin(id); C=Cin(id); fac0 = fac0(id); C_log = 10*log10(C);
y1=yin(id); y = 10*log10(y1);

%% Slope/Fraction of all data
slp = y - x; % optical/acoustic for ADV this is actually the intercept, not slope.
slpa = C_log - x; % concentration/acoustic
slpo = C./y1; % concentration/optical
fac = 100-((C-fac0)./C*100); 

xyfac = fit(slp,fac,'poly1'); 
xC = fit(fac,slpa,'power2'); 
yC = fit(fac,slpo,'exp1'); 

%% Slope/Fraction of Part 1
slp1 = y(1:29) - x(1:29); % optical/acoustic for ADV this is actually the intercept, not slope.
slpa1 = C_log(1:29) - x(1:29); % concentration/acoustic
slpo1 = C(1:29)./y1(1:29); % concentration/optical
fac1 = 100-((C(1:29)-fac0(1:29))./C(1:29)*100); 

xyfac1 = fit(slp1, fac1,'poly1'); %main.xyfac = xyfac;
xC1 = fit(fac1,slpa1,'power2'); 
yC1 = fit(fac1,slpo1,'exp1'); 

%% Slope/Fraction of Part 2
slp2 = y(30:59) - x(30:59); % optical/acoustic for ADV this is actually the intercept, not slope.
slpa2 = C_log(30:59) - x(30:59); % concentration/acoustic
slpo2 = C(30:59)./y1(30:59); % concentration/optical
fac2 = 100-((C(30:59)-fac0(30:59))./C(30:59)*100); 

xyfac2 = fit(slp2,fac2,'poly1'); 
xC2 = fit(fac2,slpa2,'power2'); 
yC2 = fit(fac2,slpo2,'exp1'); 

n=0;
for i = 1:numel(slp)
    %% use all data
    Be(i) = xyfac(slp(i));
    if Be(i) > 100
        Be(i) = 100;
    end    
    if Be(i) < 0
        Be(i) = 0.0001;
    end
    Ca_temp(i) = x(i) + xC(Be(i));
    Ca(i) = 10^(Ca_temp(i)/10);
    Co(i) = yC(Be(i))*y1(i);
    
    %% use data from Sand 100
    if i>= 1 && i <= 29
    Be1(i) = xyfac1(slp1(i));
    if Be1(i) > 100
        Be1(i) = 100;
    end    
    if Be1(i) < 0
        Be1(i) = 0.0001;
    end
    Ca1_temp(i) = x(i) + xC1(Be1(i));
    Ca1(i) = 10^(Ca1_temp(i)/10);
    Co1(i) = yC1(Be1(i))*y1(i);
    end

    %% use data from Sand 200
    if i >= 30 && i <=59
    n=n+1;   
    Be2(n) = xyfac2(slp2(n));
    if  Be2(n) > 100
        Be2(n) = 100;
    end    
    if Be2(n) < 0
        Be2(n) = 0.0001;
    end
    Ca2_temp(n) = x(i) + xC2(Be2(n));
    Ca2(n) = 10^(Ca2_temp(n)/10);
    Co2(n) = yC2(Be2(n))*y1(i);
    end
end

%% use fitted function to calculate S1S2
for j = 1:6
    slp_cal(j) = y_cal_log(j) - x_cal(j);
    Be_cal(j) = xyfac(slp_cal(j));
    if Be_cal(j) > 100
        Be_cal(j) = 100;
    end    
    if Be_cal(j) < 0
        Be_cal(j) = 0.0001;
    end
    Ca_cal_temp(j) = x_cal(j) + xC(Be_cal(j));
    Ca_cal(j) = 10^(Ca_cal_temp(j)/10);
    Co_cal(j) = yC(Be_cal(j))*y_cal(j);  
%%%%%%%%%%%%%%%%%%%%%%
    Be_cal1(j) = xyfac1(slp_cal(j));
    if  Be_cal1(j) > 100
        Be_cal1(j) = 100;
    end    
    if  Be_cal1(j) < 0
        Be_cal1(j) = 0.0001;
    end
    Ca_cal1_temp(j) = x_cal(j) + xC1(Be_cal1(j));
    Ca_cal1(j) = 10^(Ca_cal1_temp(j)/10);
    Co_cal1(j) = yC1(Be_cal1(j))*y_cal(j);  
%%%%%%%%%%%%%%%%%%%%%%
    Be_cal2(j) = xyfac2(slp_cal(j));
    if Be_cal2(j) > 100
        Be_cal2(j) = 100;
    end    
    if Be_cal2(j) < 0
        Be_cal2(j) = 0.0001;
    end
    Ca_cal2_temp(j) = x_cal(j) + xC2(Be_cal2(j));
    Ca_cal2(j) = 10^(Ca_cal2_temp(j)/10);
    Co_cal2(j) = yC2(Be_cal2(j))*y_cal(j);      
end

%% save results
main.slp = slp; main.slpa = slpa; main.slpo = slpo; 
main.xyfac = xyfac; main.xC = xC; main.yC = yC; 
main.Co_xyfac = coeffvalues(xyfac); main.Co_xC = coeffvalues(xC); 
main.Co_yC = coeffvalues(yC); 

main.slp1 = slp1; main.slpa1 = slpa1; main.slpo1 = slpo1; 
main.xyfac1 = xyfac1; main.xC1 = xC1; main.yC1 = yC1;
main.Co_xyfac1 = coeffvalues(xyfac1); main.Co_xC1 = coeffvalues(xC1); 
main.Co_yC1 = coeffvalues(yC1);

main.slp2 = slp2; main.slpa2 = slpa2; main.slpo2 = slpo2; 
main.xyfac2 = xyfac2; main.xC2 = xC2; main.yC2 = yC2;
main.Co_xyfac2 = coeffvalues(xyfac2); main.Co_xC2 = coeffvalues(xC2); 
main.Co_yC2 = coeffvalues(yC2);

main.C = C;
main.Ca = Ca'; main.Co = Co';
temp1 = [Ca1,Ca2]; temp2 = [Co1,Co2];
main.Ca12 = temp1'; main.Co12 = temp2';

main.Be = Be';
main.diff_all = rmsebias(main.C,main.Ca,main.Co);
main.diff_S1 = rmsebias(main.C(1:29),main.Ca12(1:29),main.Co12(1:29));
main.diff_S2 = rmsebias(main.C(30:59),main.Ca12(30:59),main.Co12(30:59));

main.v.Be_cal = Be_cal'; 
main.v.Ca_cal = Ca_cal'; 
main.v.Co_cal = Co_cal'; 

main.v.Be_cal1 = Be_cal1'; 
main.v.Ca_cal1 = Ca_cal1'; 
main.v.Co_cal1 = Co_cal1'; 

main.v.Be_cal2 = Be_cal2';
main.v.Ca_cal2 = Ca_cal2';
main.v.Co_cal2 = Co_cal2';

main.v.diff_all = rmsebias(data.CBev,main.v.Ca_cal, main.v.Co_cal);
main.v.diff_S1  = rmsebias(data.CBev,main.v.Ca_cal1,main.v.Co_cal1);
main.v.diff_S2  = rmsebias(data.CBev,main.v.Ca_cal2,main.v.Co_cal2);

main.diff_all.S1 = rmsebias(main.C(1:29),main.Ca(1:29),main.Co(1:29));
main.diff_all.S2 = rmsebias(main.C(30:59),main.Ca(30:59),main.Co(30:59));
main.diff_all.Be1 = rmsebias(main.Be(1:29),fac(1:29),fac(1:29));
main.diff_all.Be2 = rmsebias(main.Be(30:59),fac(30:59),fac(30:59));

figure
subplot(1,2,1)
plot(main.diff_all.Ca,'o')
title('Ca % -- ADV Hydro4')

subplot(1,2,2)
plot(main.diff_all.Co,'s')
title('Co % -- ADV Hydro4')


% save(name,'main')

% % % xyspace = linspace(min(slp), max(slp), 150);  
% % % facspace = linspace(min(fac), max(fac), 150); 
% % % 
% % % xyspace2 = linspace(min(slp), max(slp), 150);  
% % % facspace2 = linspace(min(fac2), max(fac2), 150); 
% % % facspace1 = linspace(min(fac1), max(fac1), 150); 
% % % 
% % % figure
% % % subplot(3,2,1)
% % % plot(xin(1:30),10*log10(yin(1:30)),'o','MarkerSize',8.5,'MarkerFaceColor','b','MarkerEdgeColor','b')
% % % hold on
% % % plot(xin(31:60),10*log10(yin(1:30)),'s','MarkerSize',8.5,'MarkerFaceColor','r','MarkerEdgeColor','r')
% % % title('Acoustic vs. Optical')
% % % xlabel([xlab '_{signal}'])
% % % ylabel(['10log10(' ylab '_{signal})'])
% % % grid on
% % % 
% % % subplot(3,2,2)
% % % plot(slp(1:29),fac(1:29),'o','MarkerSize',8.5,'MarkerFaceColor','b','MarkerEdgeColor','b')
% % % hold on
% % % plot(slp(30:59),fac(30:59),'s','MarkerSize',8.5,'MarkerFaceColor','r','MarkerEdgeColor','r')
% % % hold on
% % % plot(xyspace,xyfac(xyspace),'LineWidth',2.3,'color',[0 0 0])
% % % hold on
% % % plot(xyspace,xyfac1(xyspace),'LineWidth',1.7,'color','b')
% % % hold on
% % % plot(xyspace,xyfac2(xyspace),'LineWidth',1.7,'color','r')
% % % title(['From ' xlab '_{signal} / ' ylab '_{signal} => Be fraction' ])
% % % xlabel(['10log10(' ylab '_{signal}) - ' xlab '_{signal}'])
% % % ylabel('Be fraction')
% % % ylim([0 100])
% % % grid on
% % % 
% % % subplot(3,2,3)
% % % plot(xin(1:30),10*log10(Cin(1:30)),'o','MarkerSize',8.5,'MarkerFaceColor','b','MarkerEdgeColor','b')
% % % hold on
% % % plot(xin(31:60),10*log10(Cin(31:60)),'s','MarkerSize',8.5,'MarkerFaceColor','r','MarkerEdgeColor','r')
% % % title('Acoustic vs. C')
% % % xlabel([xlab '_{signal}'])
% % % ylabel('10log10(C)')
% % % grid on
% % % 
% % % subplot(3,2,4)
% % % plot(fac(1:29),slpa(1:29),'o','MarkerSize',8.5,'MarkerFaceColor','b','MarkerEdgeColor','b')
% % % hold on
% % % plot(fac(30:59),slpa(30:59),'s','MarkerSize',8.5,'MarkerFaceColor','r','MarkerEdgeColor','r')
% % % hold on
% % % plot(facspace,xC(facspace),'LineWidth',2.3,'color',[0 0 0])
% % % hold on
% % % plot(facspace,xC1(facspace),'LineWidth',1.7,'color','b')
% % % hold on
% % % plot(facspace,xC2(facspace),'LineWidth',1.7,'color','r')
% % % title(['From Be fraction => 10log10(C) - ' xlab '_{signal}, then C'])
% % % xlabel('Be fraction')
% % % ylabel(['10log10(C) - ' xlab '_{signal}'])
% % % xlim([0 100])
% % % grid on
% % % 
% % % subplot(3,2,5)
% % % plot(yin(1:30),Cin(1:30),'o','MarkerSize',8.5,'MarkerFaceColor','b','MarkerEdgeColor','b')
% % % hold on
% % % plot(yin(31:60),Cin(31:60),'s','MarkerSize',8.5,'MarkerFaceColor','r','MarkerEdgeColor','r')
% % % title('Optical vs. C')
% % % xlabel([ylab '_{signal}'])
% % % ylabel('C')
% % % grid on
% % % 
% % % subplot(3,2,6)
% % % plot(fac(1:29),slpo(1:29),'o','MarkerSize',8.5,'MarkerFaceColor','b','MarkerEdgeColor','b')
% % % hold on
% % % plot(fac(30:59),slpo(30:59),'s','MarkerSize',8.5,'MarkerFaceColor','r','MarkerEdgeColor','r')
% % % hold on
% % % plot(facspace,yC(facspace),'LineWidth',2.3,'color',[0 0 0])
% % % hold on
% % % plot(facspace,yC1(facspace),'LineWidth',1.7,'color','b')
% % % hold on
% % % plot(facspace,yC2(facspace),'LineWidth',1.7,'color','r')
% % % title(['From Be fraction => C/' ylab '_{signal}, then C'])
% % % xlabel('Be fraction')
% % % ylabel(['C / ' ylab '_{signal}'])
% % % xlim([0 100])
% % % grid on