%
% Plotting script for dynamic quantites from dynasor output file.
%

clear all
close all
clc

% Colors
red=[0.9 0.3 0.3];
blue=[0.3 0.3 0.9];
black=[0.1 0.1 0.1];
green=[0.3 0.9 0.3];
purple=[0.8 0.4 0.8];
% Colors

fontsize=20;
fontstyle='Times';

% Conversion factors
hbar=1.05457173*1e-34; % m2 kg / s
J2ev=6.24150934e18;

% Load dynasor output
dynasor_outT1400_dynamical_old

w=hbar*w*J2ev*1000*1e15; %meV

%% F(k,t) S(k,w) C(k,t) C(k,w)

kind1=7; % Two k-values
kind2=22;

wMax=80;
tMax=900;
LW=5.0;


figure('Position',[100 20,1500,800],'Color','w')

subplot(2,2,1)
plot(t,F_k_t_0_0(:,kind2),'Color',red,'LineWidth',LW)
hold on
plot(t,F_k_t_0_0(:,kind1),'Color',blue,'LineWidth',LW)
ylabel('F(k,t)')
h=legend(strcat('k=',num2str(k(kind2)),' nm^{-1}'),strcat('k=',num2str(k(kind1)),' nm^{-1}'));
set(h,'Edgecolor','w')
xlim([0 tMax])
ylim([-5 20]*1e-3)

subplot(2,2,2)
plot(w(2:end),S_k_w_0_0(2:end,kind2),'Color',red,'LineWidth',LW)
hold on
plot(w(2:end),S_k_w_0_0(2:end,kind1),'Color',blue,'LineWidth',LW)
h=legend(strcat('k=',num2str(k(kind2)),' nm^{-1}'),strcat('k=',num2str(k(kind1)),' nm^{-1}'));
set(h,'Edgecolor','w')
ylabel('S(k,w)')
xlim([0 wMax])
ylim([0 8])

subplot(2,2,3)
plot(t,Cl_k_t_0_0(:,kind2),'Color',red,'LineWidth',LW)
hold on
plot(t,Ct_k_t_0_0(:,kind2),'Color',blue,'LineWidth',LW)
xlabel('Time [fs]')
ylabel('C(k,t)')
h=legend(strcat('Longitudinal k=',num2str(k(kind2)),' nm^{-1}'),strcat('Transversal k=',num2str(k(kind2)),' nm^{-1}'));
set(h,'Edgecolor','w')
xlim([0 tMax])
ylim([-0.3 0.15])

subplot(2,2,4)
plot(w,Cl_k_w_0_0(:,kind2),'Color',red,'LineWidth',LW)
hold on
plot(w,Ct_k_w_0_0(:,kind2),'Color',blue,'LineWidth',LW)
h=legend(strcat('Longitudinal k=',num2str(k(kind2)),' nm^{-1}'),strcat('Transversal k=',num2str(k(kind2)),' nm^{-1}'));
set(h,'Edgecolor','w')
xlabel('\omega [meV]')
ylabel('C(k,w)')
xlim([0 wMax])
ylim([0 70])


set(findall(gcf,'-property','FontSize'),'FontSize',fontsize)
set(findall(gcf,'-property','FontName'),'FontName',fontstyle)

%% map of Cl(k,w) and Ct(k,w)

kstart=8;
wMax=80;
[tmp,wInd]=min(abs(w-wMax));
xlims=[k(kstart) k(end)];


figure('Position',[100 400,1500,650],'Color','w')

subplot(1,2,1)
surf(k(kstart:end),w(1:wInd),Cl_k_w_0_0(1:wInd,kstart:end),'EdgeColor','none','LineStyle','none')
xlim(xlims)
ylim([0 w(wInd)])
ylabel('\omega [meV]')
xlabel('k [nm^{-1}] ')
title('Cl(k,\omega)')
view([0 90])

subplot(1,2,2)
surf(k(kstart:end),w(1:wInd),Ct_k_w_0_0(1:wInd,kstart:end),'EdgeColor','none','LineStyle','none')
xlim(xlims)
ylim([0 w(wInd)])
xlabel('k [nm^{-1}] ')
title('Ct(k,\omega)')
view([0 90])

set(findall(gcf,'-property','FontSize'),'FontSize',fontsize)
set(findall(gcf,'-property','FontName'),'FontName',fontstyle)

%% smooth map of S(k,w) and Cl(k,w)

kstart=8;
wMax=80;
[tmp,wInd]=min(abs(w-wMax));
xlims=[k(kstart) k(end)];

ave=15;
smooth_CL=zeros(wInd,length(k)-kstart+1);
smooth_CT=zeros(wInd,length(k)-kstart+1);
for i=kstart:length(k)
    smooth_CL(:,i)=smooth(Cl_k_w_0_0(1:wInd,i),ave,'lowess');
    smooth_CT(:,i)=smooth(Ct_k_w_0_0(1:wInd,i),ave,'lowess');
end


figure('Position',[100 400,1500,650],'Color','w')

subplot(1,2,1)
surf(k(kstart:end),w(1:wInd),smooth_CL(1:wInd,kstart:end),'EdgeColor','none','LineStyle','none')
xlim(xlims)
ylim([0 w(wInd)])
ylabel('\omega [meV]')
xlabel('k [nm^{-1}] ')
title('smoothed  Cl(k,\omega)')
view([0 90])

subplot(1,2,2)
surf(k(kstart:end),w(1:wInd),smooth_CT(1:wInd,kstart:end),'EdgeColor','none','LineStyle','none')
xlim(xlims)
ylim([0 w(wInd)])
xlabel('k [nm^{-1}] ')
title('smoothed  Ct(k,\omega)')
view([0 90])

set(findall(gcf,'-property','FontSize'),'FontSize',fontsize)
set(findall(gcf,'-property','FontName'),'FontName',fontstyle)
