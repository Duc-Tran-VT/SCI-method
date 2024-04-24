% function main = facAQUA(xin,yin,Cin)
clear all; clc; %close all;

load('All.mat'); load('S1S2.mat');

xin = all.aqua4.^2; x_cal = data.aqua4.^2; xlab = 'Hydro4';
yin = all.hydro4; y_cal = data.hydro4; ylab = 'Hydro4'; name = 'AQUA4_Hydro4';
Cin = all.C;

fac0 = 0.9*[15,25,50,100,150,200,11.25,18.75,37.5,75,112.5,150,7.5,12.5,25,...
    50,75,100,3.75,6.25,12.5,25,37.5,50,0.001,0.001,0.001,0.001,0.001,0.001,...
        15,25,50,100,150,200,11.25,18.75,37.5,75,112.5,150,7.5,12.5,25,...
    50,75,100,3.75,6.25,12.5,25,37.5,50,0.001,0.001,0.001,0.001,0.001,0.001]'; 

id = [1:11 13:45 49:60]'; % part 1 from 1:29, part 2 from 30:56 (in new index)
x=xin(id); C=Cin(id); fac0 = fac0(id); y=yin(id);

%% Slope/Fraction of all data
slp = y./x; % optical/acoustic for ADV this is actually the intercept, not slope.
slpa = C./x; % concentration/acoustic
slpo = C./y; % concentration/optical
fac = 100-((C-fac0)./C*100);

xyfac = fit(slp,fac,'power2');
xC = fit(fac,slpa,'power2');
yC = fit(fac,slpo,'exp1');

%% Slope/Fraction of Part 1
slp1 = y(1:29)./x(1:29); % optical/acoustic for ADV this is actually the intercept, not slope.
slpa1 = C(1:29)./x(1:29); % concentration/acoustic
slpo1 = C(1:29)./y(1:29); % concentration/optical
fac1 = 100-((C(1:29)-fac0(1:29))./C(1:29)*100);

xyfac1 = fit(slp1, fac1, 'power2'); %main.xyfac = xyfac;
xC1 = fit(fac1,slpa1,'power2');
yC1 = fit(fac1,slpo1,'exp1');

%% Slope/Fraction of Part 2
slp2 = y(30:56)./x(30:56); % optical/acoustic for ADV this is actually the intercept, not slope.
slpa2 = C(30:56)./x(30:56); % concentration/acoustic
slpo2 = C(30:56)./y(30:56); % concentration/optical
fac2 = 100-((C(30:56)-fac0(30:56))./C(30:56)*100);

xyfac2 = fit(slp2,fac2,'power2');
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
    Ca(i) = xC(Be(i))*x(i);
    Co(i) = yC(Be(i))*y(i);
    
    %% use data from Sand 100
    if i>= 1 && i <= 29
        Be1(i) = xyfac1(slp1(i));
        if Be1(i) > 100
            Be1(i) = 100;
        end
        if Be1(i) < 0
            Be1(i) = 0.0001;
        end
        Ca1(i) = xC1(Be1(i))*x(i);
        Co1(i) = yC1(Be1(i))*y(i);
    end
    
    %% use data from Sand 200
    if i >= 30 && i <=56
        n=n+1;
        Be2(n) = xyfac2(slp2(n));
        if Be2(n) > 100
            Be2(n) = 100;
        end
        if Be2(n) < 0
            Be2(n) = 0.0001;
        end
        Ca2(n) = xC2(Be2(n))*x(i);
        Co2(n) = yC2(Be2(n))*y(i);
    end
end

%% use fitted function to calculate S1S2
for j = 1:6
    slp_cal(j) = y_cal(j)./x_cal(j);
    Be_cal(j) = xyfac(slp_cal(j));
    if Be_cal(j) > 100
        Be_cal(j) = 100;
    end
    if Be_cal(j) < 0
        Be_cal(j) = 0.0001;
    end
    Ca_cal(j) = xC(Be_cal(j))*x_cal(j);
    Co_cal(j) = yC(Be_cal(j))*y_cal(j);
    %%%%%%%%%%%%%%%%%%%%%%
    Be_cal1(j) = xyfac1(slp_cal(j));
    if Be_cal1(j) > 100
        Be_cal1(j) = 100;
    end
    if Be_cal1(j) < 0
        Be_cal1(j) = 0.0001;
    end
    Ca_cal1(j) = xC1(Be_cal1(j))*x_cal(j);
    Co_cal1(j) = yC1(Be_cal1(j))*y_cal(j);
    %%%%%%%%%%%%%%%%%%%%%%
    Be_cal2(j) = xyfac2(slp_cal(j));
    if Be_cal2(j) > 100
        Be_cal2(j) = 100;
    end
    if Be_cal2(j) < 0
        Be_cal2(j) = 0.0001;
    end
    Ca_cal2(j) = xC2(Be_cal2(j))*x_cal(j);
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
% 
% main.C = C;
% main.Ca = Ca'; main.Co = Co';
% temp1 = [Ca1,Ca2]; temp2 = [Co1,Co2];
% main.Ca12 = temp1'; main.Co12 = temp2';
% 
% main.Be = Be';
% 
% main.diff_all = rmsebias(main.C,main.Ca,main.Co);
% main.diff_S1  = rmsebias(main.C(1:29),main.Ca12(1:29),main.Co12(1:29));
% main.diff_S2  = rmsebias(main.C(30:56),main.Ca12(30:56),main.Co12(30:56));
% 
% main.v.Be_cal = Be_cal'; 
% main.v.Ca_cal = Ca_cal'; 
% main.v.Co_cal = Co_cal'; 
% 
% main.v.Be_cal1 = Be_cal1'; 
% main.v.Ca_cal1 = Ca_cal1'; 
% main.v.Co_cal1 = Co_cal1'; 
% 
% main.v.Be_cal2 = Be_cal2';
% main.v.Ca_cal2 = Ca_cal2';
% main.v.Co_cal2 = Co_cal2';
% 
% main.v.diff_all = rmsebias(data.CBev,main.v.Ca_cal, main.v.Co_cal);
% main.v.diff_S1  = rmsebias(data.CBev,main.v.Ca_cal1,main.v.Co_cal1);
% main.v.diff_S2  = rmsebias(data.CBev,main.v.Ca_cal2,main.v.Co_cal2);
% 
% main.diff_all.S1 = rmsebias(main.C(1:29),main.Ca(1:29),main.Co(1:29));
% main.diff_all.S2 = rmsebias(main.C(30:56),main.Ca(30:56),main.Co(30:56));
% main.diff_all.Be1 = rmsebias(main.Be(1:29),fac(1:29),fac(1:29));
% main.diff_all.Be2 = rmsebias(main.Be(30:56),fac(30:56),fac(30:56));
% 
% figure
% subplot(1,2,1)
% plot(main.diff_all.Ca,'o')
% title('Ca % -- AQUA3 Hydro4')
% 
% subplot(1,2,2)
% plot(main.diff_all.Co,'s')
% title('Co % -- AQUA3 Hydro4')


% save(name,'main')

%
%% Figure Sand 1
xyspace = linspace(min(slp), max(slp), 150);  
facspace = linspace(min(fac), max(fac), 150); 

xyspace2 = linspace(min(slp), max(slp), 150);  
facspace2 = linspace(min(fac2), max(fac2), 150); 
facspace1 = linspace(min(fac1), max(fac1), 150); 

figure
subplot(3,2,1)
plot(xin(1:6),yin(1:6),':o','MarkerSize',8.5)
hold on
plot(xin(7:11),yin(7:11),':o','MarkerSize',8.5,'MarkerFaceColor','b','MarkerEdgeColor','b')
hold on
plot(xin(12:17),yin(12:17),':o','MarkerSize',8.5,'MarkerFaceColor','g','MarkerEdgeColor','g')
hold on
plot(xin(18:23),yin(18:23),':o','MarkerSize',8.5,'MarkerFaceColor','r','MarkerEdgeColor','r')
hold on
plot(xin(24:29),yin(24:29),':o','MarkerSize',8.5,'MarkerFaceColor','k','MarkerEdgeColor','k')
hold on
plot(xin(30:35),yin(30:35),':s','MarkerSize',8.5)
hold on
plot(xin(36:41),yin(36:41),':s','MarkerSize',8.5,'MarkerFaceColor','b','MarkerEdgeColor','g')
hold on
plot(xin(42:44),yin(42:44),':s','MarkerSize',8.5,'MarkerFaceColor','g','MarkerEdgeColor','r')
hold on
plot(xin(45:50),yin(45:50),':s','MarkerSize',8.5,'MarkerFaceColor','r','MarkerEdgeColor','k')
hold on
plot(xin(51:56),yin(51:56),':s','MarkerSize',8.5,'MarkerFaceColor','k','MarkerEdgeColor','k')
hold on
plot(x_cal(1),y_cal(1),'v','MarkerSize',8.5)
hold on
plot(x_cal(2),y_cal(2),'v','MarkerSize',8.5,'MarkerFaceColor','b','MarkerEdgeColor','b')
hold on
plot(x_cal(3),y_cal(3),'v','MarkerSize',8.5,'MarkerFaceColor','g','MarkerEdgeColor','g')
hold on
plot(x_cal(4),y_cal(4),'v','MarkerSize',8.5,'MarkerFaceColor','r','MarkerEdgeColor','r')
hold on
plot(x_cal(5),y_cal(5),'v','MarkerSize',8.5,'MarkerFaceColor','k','MarkerEdgeColor','k')
hold on
plot(x_cal(6),y_cal(6),'v','MarkerSize',12,'MarkerFaceColor','b','MarkerEdgeColor','b')
title('Acoustic vs. Optical', 'FontSize',15)
xlabel(xlab, 'FontSize',15)
ylabel(ylab, 'FontSize',15)
text(10,16,'a)','FontSize',22)
grid on

subplot(3,2,2)
plot(slp(1:6),fac(1:6),'o','MarkerSize',8.5)
hold on
plot(slp(7:11),fac(7:11),'o','MarkerSize',8.5,'MarkerFaceColor','b','MarkerEdgeColor','b')
hold on
plot(slp(12:17),fac(12:17),'o','MarkerSize',8.5,'MarkerFaceColor','g','MarkerEdgeColor','g')
hold on
plot(slp(18:23),fac(18:23),'o','MarkerSize',8.5,'MarkerFaceColor','r','MarkerEdgeColor','r')
hold on
plot(slp(24:29),fac(24:29),'o','MarkerSize',8.5,'MarkerFaceColor','k','MarkerEdgeColor','k')
hold on
plot(slp(30:35),fac(30:35),'s','MarkerSize',8.5)
hold on
plot(slp(36:41),fac(36:41),'s','MarkerSize',8.5,'MarkerFaceColor','b','MarkerEdgeColor','g')
hold on
plot(slp(42:44),fac(42:44),'s','MarkerSize',8.5,'MarkerFaceColor','g','MarkerEdgeColor','r')
hold on
plot(slp(45:50),fac(45:50),'s','MarkerSize',8.5,'MarkerFaceColor','r','MarkerEdgeColor','k')
hold on
plot(slp(51:56),fac(51:56),'s','MarkerSize',8.5,'MarkerFaceColor','k','MarkerEdgeColor','k')
hold on
plot(y_cal(1)/x_cal(1),data.Bev(1),'v','MarkerSize',8.5)
hold on
plot(y_cal(2)/x_cal(2),data.Bev(2),'v','MarkerSize',8.5,'MarkerFaceColor','b','MarkerEdgeColor','b')
hold on
plot(y_cal(3)/x_cal(3),data.Bev(3),'v','MarkerSize',8.5,'MarkerFaceColor','g','MarkerEdgeColor','g')
hold on
plot(y_cal(4)/x_cal(4),data.Bev(4),'v','MarkerSize',8.5,'MarkerFaceColor','r','MarkerEdgeColor','r')
hold on
plot(y_cal(5)/x_cal(5),data.Bev(5),'v','MarkerSize',8.5,'MarkerFaceColor','k','MarkerEdgeColor','k')
hold on
plot(y_cal(6)/x_cal(6),data.Bev(6),'v','MarkerSize',12,'MarkerFaceColor','b','MarkerEdgeColor','b')
plot(xyspace,xyfac(xyspace),'LineWidth',1.7,'color','b')
legend('Be','BeS1_{31}','BeS1_{11}','BeS1_{13}','S1','fitted line','FontSize',15)
title(['Step1: Be_{fraction} = f(' xlab ' / ' ylab ')' ], 'FontSize',15)
xlabel([ylab ' / ' xlab ], 'FontSize',15)
ylabel('Be_{fraction} [%]', 'FontSize',15)
ylim([0 100])
text(0.1,85,'b)','FontSize',22)
grid on

subplot(3,2,3)
plot(xin(1:6),Cin(1:6),':o','MarkerSize',8.5)
hold on
plot(xin(7:11),Cin(7:11),':o','MarkerSize',8.5,'MarkerFaceColor','b','MarkerEdgeColor','b')
hold on
plot(xin(12:17),Cin(12:17),':o','MarkerSize',8.5,'MarkerFaceColor','g','MarkerEdgeColor','g')
hold on
plot(xin(18:23),Cin(18:23),':o','MarkerSize',8.5,'MarkerFaceColor','r','MarkerEdgeColor','r')
hold on
plot(xin(24:29),Cin(24:29),':o','MarkerSize',8.5,'MarkerFaceColor','k','MarkerEdgeColor','k')
hold on
plot(xin(30:35),Cin(30:35),':s','MarkerSize',8.5)
hold on
plot(xin(36:41),Cin(36:41),':s','MarkerSize',8.5,'MarkerFaceColor','b','MarkerEdgeColor','g')
hold on
plot(xin(42:44),Cin(42:44),':s','MarkerSize',8.5,'MarkerFaceColor','g','MarkerEdgeColor','r')
hold on
plot(xin(45:50),Cin(45:50),':s','MarkerSize',8.5,'MarkerFaceColor','r','MarkerEdgeColor','k')
hold on
plot(xin(51:56),Cin(51:56),':s','MarkerSize',8.5,'MarkerFaceColor','k','MarkerEdgeColor','k')
title('Acoustic vs. C', 'FontSize',15)
xlabel(xlab, 'FontSize',15)
ylabel('C [mg/L]', 'FontSize',15)
text(9,160,'c)','FontSize',22)
grid on

subplot(3,2,4)
plot(fac(1:6),slpa(1:6),'o','MarkerSize',8.5)
hold on
plot(fac(7:11),slpa(7:11),'o','MarkerSize',8.5,'MarkerFaceColor','b','MarkerEdgeColor','b')
hold on
plot(fac(12:17),slpa(12:17),'o','MarkerSize',8.5,'MarkerFaceColor','g','MarkerEdgeColor','g')
hold on
plot(fac(18:23),slpa(18:23),'o','MarkerSize',8.5,'MarkerFaceColor','r','MarkerEdgeColor','r')
hold on
plot(fac(24:29),slpa(24:29),'o','MarkerSize',8.5,'MarkerFaceColor','k','MarkerEdgeColor','k')
hold on
plot(fac(30:35),slpa(30:35),'s','MarkerSize',8.5)
hold on
plot(fac(36:41),slpa(36:41),'s','MarkerSize',8.5,'MarkerFaceColor','b','MarkerEdgeColor','g')
hold on
plot(fac(42:44),slpa(42:44),'s','MarkerSize',8.5,'MarkerFaceColor','g','MarkerEdgeColor','r')
hold on
plot(fac(45:50),slpa(45:50),'s','MarkerSize',8.5,'MarkerFaceColor','r','MarkerEdgeColor','k')
hold on
plot(fac(51:56),slpa(51:56),'s','MarkerSize',8.5,'MarkerFaceColor','k','MarkerEdgeColor','k')
hold on
plot(facspace,xC(facspace),'LineWidth',1.7,'color','b')
title(['Step 2a: C / ' xlab ' = f(Be_{fraction})'], 'FontSize',15)
xlabel('Be_{fraction} [%]', 'FontSize',15)
ylabel(['C / ' xlab], 'FontSize',15)
xlim([0 100])
text(10,5.3,'d)','FontSize',22)
grid on

subplot(3,2,5)
plot(yin(1:6),Cin(1:6),':o','MarkerSize',8.5)
hold on
plot(yin(7:12),Cin(7:12),':o','MarkerSize',8.5,'MarkerFaceColor','b','MarkerEdgeColor','b')
hold on
plot(yin(13:18),Cin(13:18),':o','MarkerSize',8.5,'MarkerFaceColor','g','MarkerEdgeColor','g')
hold on
plot(yin(19:24),Cin(19:24),':o','MarkerSize',8.5,'MarkerFaceColor','r','MarkerEdgeColor','r')
hold on
plot(yin(25:30),Cin(25:30),':o','MarkerSize',8.5,'MarkerFaceColor','k','MarkerEdgeColor','k')
title('Optical vs. C', 'FontSize',15)
xlabel(ylab, 'FontSize',15)
ylabel('C [mg/L]', 'FontSize',15)
text(1.5,160,'e)','FontSize',22)
grid on

subplot(3,2,6)
plot(fac(1:6),slpo(1:6),'o','MarkerSize',8.5)
hold on
plot(fac(7:11),slpo(7:11),'o','MarkerSize',8.5,'MarkerFaceColor','b','MarkerEdgeColor','b')
hold on
plot(fac(12:17),slpo(12:17),'o','MarkerSize',8.5,'MarkerFaceColor','g','MarkerEdgeColor','g')
hold on
plot(fac(18:23),slpo(18:23),'o','MarkerSize',8.5,'MarkerFaceColor','r','MarkerEdgeColor','r')
hold on
plot(fac(24:29),slpo(24:29),'o','MarkerSize',8.5,'MarkerFaceColor','k','MarkerEdgeColor','k')
hold on
plot(fac(30:35),slpo(30:35),'s','MarkerSize',8.5)
hold on
plot(fac(36:41),slpo(36:41),'s','MarkerSize',8.5,'MarkerFaceColor','b','MarkerEdgeColor','g')
hold on
plot(fac(42:44),slpo(42:44),'s','MarkerSize',8.5,'MarkerFaceColor','g','MarkerEdgeColor','r')
hold on
plot(fac(45:50),slpo(45:50),'s','MarkerSize',8.5,'MarkerFaceColor','r','MarkerEdgeColor','k')
hold on
plot(fac(51:56),slpo(51:56),'s','MarkerSize',8.5,'MarkerFaceColor','k','MarkerEdgeColor','k')
hold on
plot(facspace,yC(facspace),'LineWidth',1.7,'color','b')
title(['Step 2b: C/ ' ylab ' = f(Be_{fraction})'], 'FontSize',15)
xlabel('Be_{fraction} [%]', 'FontSize',15)
ylabel(['C / ' ylab], 'FontSize',15)
xlim([0 100])
text(10,26,'f)','FontSize',22)
grid on

%%%%%%%%%%%%%
figure
subplot(3,2,1)
plot(xin(31:36),yin(31:36),':o','MarkerSize',8.5)
hold on
plot(xin(37:42),yin(37:42),':o','MarkerSize',8.5,'MarkerFaceColor','b','MarkerEdgeColor','b')
hold on
plot(xin(43:48),yin(43:48),':o','MarkerSize',8.5,'MarkerFaceColor','g','MarkerEdgeColor','g')
hold on
plot(xin(49:54),yin(49:54),':o','MarkerSize',8.5,'MarkerFaceColor','r','MarkerEdgeColor','r')
hold on
plot(xin(55:60),yin(55:60),':o','MarkerSize',8.5,'MarkerFaceColor','k','MarkerEdgeColor','k')
title('Acoustic vs. Optical', 'FontSize',15)
xlabel(xlab, 'FontSize',15)
ylabel(ylab, 'FontSize',15)
text(10,16,'a)','FontSize',22)
grid on

subplot(3,2,2)
plot(slp1(1:6),fac1(1:6),'o','MarkerSize',8.5)
hold on
plot(slp1(7:11),fac1(7:11),'o','MarkerSize',8.5,'MarkerFaceColor','b','MarkerEdgeColor','b')
hold on
plot(slp1(12:17),fac1(12:17),'o','MarkerSize',8.5,'MarkerFaceColor','g','MarkerEdgeColor','g')
hold on
plot(slp1(18:23),fac1(18:23),'o','MarkerSize',8.5,'MarkerFaceColor','r','MarkerEdgeColor','r')
hold on
plot(slp1(24:29),fac1(24:29),'o','MarkerSize',8.5,'MarkerFaceColor','k','MarkerEdgeColor','k')
hold on
plot(xyspace,xyfac1(xyspace),'LineWidth',1.7,'color','b')
legend('Be','BeS1_{31}','BeS1_{11}','BeS1_{13}','S1','fitted line','FontSize',15)
title(['Step1: Be_{fraction} = f(' xlab ' / ' ylab ')' ], 'FontSize',15)
xlabel([ylab ' / ' xlab ], 'FontSize',15)
ylabel('Be_{fraction} [%]', 'FontSize',15)
ylim([0 100])
text(0.1,85,'b)','FontSize',22)
grid on

subplot(3,2,3)
plot(xin(31:36),Cin(31:36),':o','MarkerSize',8.5)
hold on
plot(xin(37:42),Cin(37:42),':o','MarkerSize',8.5,'MarkerFaceColor','b','MarkerEdgeColor','b')
hold on
plot(xin(43:48),Cin(43:48),':o','MarkerSize',8.5,'MarkerFaceColor','g','MarkerEdgeColor','g')
hold on
plot(xin(49:54),Cin(49:54),':o','MarkerSize',8.5,'MarkerFaceColor','r','MarkerEdgeColor','r')
hold on
plot(xin(55:60),Cin(55:60),':o','MarkerSize',8.5,'MarkerFaceColor','k','MarkerEdgeColor','k')
title('Acoustic vs. C', 'FontSize',15)
xlabel(xlab, 'FontSize',15)
ylabel('C [mg/L]', 'FontSize',15)
text(9,160,'c)','FontSize',22)
grid on

subplot(3,2,4)
plot(fac1(1:6),slpa1(1:6),'o','MarkerSize',8.5)
hold on
plot(fac1(7:11),slpa1(7:11),'o','MarkerSize',8.5,'MarkerFaceColor','b','MarkerEdgeColor','b')
hold on
plot(fac1(12:17),slpa1(12:17),'o','MarkerSize',8.5,'MarkerFaceColor','g','MarkerEdgeColor','g')
hold on
plot(fac1(18:23),slpa1(18:23),'o','MarkerSize',8.5,'MarkerFaceColor','r','MarkerEdgeColor','r')
hold on
plot(fac1(24:29),slpa1(24:29),'o','MarkerSize',8.5,'MarkerFaceColor','k','MarkerEdgeColor','k')
hold on
plot(facspace,xC1(facspace),'LineWidth',1.7,'color','b')
title(['Step 2a: C / ' xlab ' = f(Be_{fraction})'], 'FontSize',15)
xlabel('Be_{fraction} [%]', 'FontSize',15)
ylabel(['C / ' xlab], 'FontSize',15)
xlim([0 100])
text(10,5.3,'d)','FontSize',22)
grid on

subplot(3,2,5)
plot(yin(31:36),Cin(31:36),':o','MarkerSize',8.5)
hold on
plot(yin(37:42),Cin(37:42),':o','MarkerSize',8.5,'MarkerFaceColor','b','MarkerEdgeColor','b')
hold on
plot(yin(43:48),Cin(43:48),':o','MarkerSize',8.5,'MarkerFaceColor','g','MarkerEdgeColor','g')
hold on
plot(yin(49:54),Cin(49:54),':o','MarkerSize',8.5,'MarkerFaceColor','r','MarkerEdgeColor','r')
hold on
plot(yin(55:60),Cin(55:60),':o','MarkerSize',8.5,'MarkerFaceColor','k','MarkerEdgeColor','k')
title('Optical vs. C', 'FontSize',15)
xlabel(ylab, 'FontSize',15)
ylabel('C [mg/L]', 'FontSize',15)
text(1.5,160,'e)','FontSize',22)
grid on

subplot(3,2,6)
plot(fac1(1:6),slpo1(1:6),'o','MarkerSize',8.5)
hold on
plot(fac1(7:11),slpo1(7:11),'o','MarkerSize',8.5,'MarkerFaceColor','b','MarkerEdgeColor','b')
hold on
plot(fac1(12:17),slpo1(12:17),'o','MarkerSize',8.5,'MarkerFaceColor','g','MarkerEdgeColor','g')
hold on
plot(fac1(18:23),slpo1(18:23),'o','MarkerSize',8.5,'MarkerFaceColor','r','MarkerEdgeColor','r')
hold on
plot(fac1(24:29),slpo1(24:29),'o','MarkerSize',8.5,'MarkerFaceColor','k','MarkerEdgeColor','k')
hold on
plot(facspace,yC1(facspace),'LineWidth',1.7,'color','b')
title(['Step 2b: C/ ' ylab ' = f(Be_{fraction})'], 'FontSize',15)
xlabel('Be_{fraction} [%]', 'FontSize',15)
ylabel(['C / ' ylab], 'FontSize',15)
xlim([0 100])
text(10,26,'f)','FontSize',22)
grid on


%%%%%%%%%%%%%
figure
subplot(3,2,1)
plot(xin(31:36),yin(31:36),':o','MarkerSize',8.5)
hold on
plot(xin(37:42),yin(37:42),':o','MarkerSize',8.5,'MarkerFaceColor','b','MarkerEdgeColor','b')
hold on
plot(xin(43:48),yin(43:48),':o','MarkerSize',8.5,'MarkerFaceColor','g','MarkerEdgeColor','g')
hold on
plot(xin(49:54),yin(49:54),':o','MarkerSize',8.5,'MarkerFaceColor','r','MarkerEdgeColor','r')
hold on
plot(xin(55:60),yin(55:60),':o','MarkerSize',8.5,'MarkerFaceColor','k','MarkerEdgeColor','k')
title('Acoustic vs. Optical', 'FontSize',15)
xlabel(xlab, 'FontSize',15)
ylabel(ylab, 'FontSize',15)
text(10,16,'a)','FontSize',22)
grid on

subplot(3,2,2)
plot(slp2(1:6),fac2(1:6),'o','MarkerSize',8.5)
hold on
plot(slp2(7:12),fac2(7:12),'o','MarkerSize',8.5,'MarkerFaceColor','b','MarkerEdgeColor','b')
hold on
plot(slp2(13:15),fac2(13:15),'o','MarkerSize',8.5,'MarkerFaceColor','g','MarkerEdgeColor','g')
hold on
plot(slp2(16:21),fac2(16:21),'o','MarkerSize',8.5,'MarkerFaceColor','r','MarkerEdgeColor','r')
hold on
plot(slp2(22:27),fac2(22:27),'o','MarkerSize',8.5,'MarkerFaceColor','k','MarkerEdgeColor','k')
hold on
plot(xyspace2,xyfac2(xyspace2),'LineWidth',1.7,'color','b')
legend('Be','BeS1_{31}','BeS1_{11}','BeS1_{13}','S1','fitted line','FontSize',15)
title(['Step1: Be_{fraction} = f(' xlab ' / ' ylab ')' ], 'FontSize',15)
xlabel([ylab ' / ' xlab ], 'FontSize',15)
ylabel('Be_{fraction} [%]', 'FontSize',15)
ylim([0 100])
text(0.1,85,'b)','FontSize',22)
grid on

subplot(3,2,3)
plot(xin(31:36),Cin(31:36),':o','MarkerSize',8.5)
hold on
plot(xin(37:42),Cin(37:42),':o','MarkerSize',8.5,'MarkerFaceColor','b','MarkerEdgeColor','b')
hold on
plot(xin(43:48),Cin(43:48),':o','MarkerSize',8.5,'MarkerFaceColor','g','MarkerEdgeColor','g')
hold on
plot(xin(49:54),Cin(49:54),':o','MarkerSize',8.5,'MarkerFaceColor','r','MarkerEdgeColor','r')
hold on
plot(xin(55:60),Cin(55:60),':o','MarkerSize',8.5,'MarkerFaceColor','k','MarkerEdgeColor','k')
title('Acoustic vs. C', 'FontSize',15)
xlabel(xlab, 'FontSize',15)
ylabel('C [mg/L]', 'FontSize',15)
text(9,160,'c)','FontSize',22)
grid on

subplot(3,2,4)
plot(fac2(1:6),slpa2(1:6),'o','MarkerSize',8.5)
hold on
plot(fac2(7:12),slpa2(7:12),'o','MarkerSize',8.5,'MarkerFaceColor','b','MarkerEdgeColor','b')
hold on
plot(fac2(13:15),slpa2(13:15),'o','MarkerSize',8.5,'MarkerFaceColor','g','MarkerEdgeColor','g')
hold on
plot(fac2(16:21),slpa2(16:21),'o','MarkerSize',8.5,'MarkerFaceColor','r','MarkerEdgeColor','r')
hold on
plot(fac2(22:27),slpa2(22:27),'o','MarkerSize',8.5,'MarkerFaceColor','k','MarkerEdgeColor','k')
hold on
plot(facspace,xC2(facspace),'LineWidth',1.7,'color','b')
title(['Step 2a: C / ' xlab ' = f(Be_{fraction})'], 'FontSize',15)
xlabel('Be_{fraction} [%]', 'FontSize',15)
ylabel(['C / ' xlab], 'FontSize',15)
xlim([0 100])
text(10,5.3,'d)','FontSize',22)
grid on

subplot(3,2,5)
plot(yin(31:36),Cin(31:36),':o','MarkerSize',8.5)
hold on
plot(yin(37:42),Cin(37:42),':o','MarkerSize',8.5,'MarkerFaceColor','b','MarkerEdgeColor','b')
hold on
plot(yin(43:48),Cin(43:48),':o','MarkerSize',8.5,'MarkerFaceColor','g','MarkerEdgeColor','g')
hold on
plot(yin(49:54),Cin(49:54),':o','MarkerSize',8.5,'MarkerFaceColor','r','MarkerEdgeColor','r')
hold on
plot(yin(55:60),Cin(55:60),':o','MarkerSize',8.5,'MarkerFaceColor','k','MarkerEdgeColor','k')
title('Optical vs. C', 'FontSize',15)
xlabel(ylab, 'FontSize',15)
ylabel('C [mg/L]', 'FontSize',15)
text(1.5,160,'e)','FontSize',22)
grid on

subplot(3,2,6)
plot(fac2(1:6),slpo2(1:6),'o','MarkerSize',8.5)
hold on
plot(fac2(7:12),slpo2(7:12),'o','MarkerSize',8.5,'MarkerFaceColor','b','MarkerEdgeColor','b')
hold on
plot(fac2(13:15),slpo2(13:15),'o','MarkerSize',8.5,'MarkerFaceColor','g','MarkerEdgeColor','g')
hold on
plot(fac2(16:21),slpo2(16:21),'o','MarkerSize',8.5,'MarkerFaceColor','r','MarkerEdgeColor','r')
hold on
plot(fac2(22:27),slpo2(22:27),'o','MarkerSize',8.5,'MarkerFaceColor','k','MarkerEdgeColor','k')
hold on
plot(facspace,yC2(facspace),'LineWidth',1.7,'color','b')
title(['Step 2b: C/ ' ylab ' = f(Be_{fraction})'], 'FontSize',15)
xlabel('Be_{fraction} [%]', 'FontSize',15)
ylabel(['C / ' ylab], 'FontSize',15)
xlim([0 100])
text(10,26,'f)','FontSize',22)
grid on