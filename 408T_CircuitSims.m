clear all
close all

%% calculate F_cyc
B = [1 1.25 1.5];
q = 1.6E-19;
m = 1.67E-27;
f_cyc = (q.*B)/(4.*pi.*m);

%% Calculate corresponding variables for resonance - Constant L 
res_freq = f_cyc;
res_omega = res_freq * (2*pi);
% C = 71E-12;
L = 2.7E-6;
R = 0.8;
% L = 1/(res_omega^2*C);
C = 1./(res_omega.^2.*L);
% res_omega = 1/sqrt(L*C)
% res_freq = res_omega/(2*pi)

freq = [5E6:1.5E2:15E6]'; %sweep from 0 to 20MHz
omega = freq*2*pi;

Z1 = R + i.*omega.*L;
Z2 = 1./(i.*omega*C);
Z = (Z1.*Z2)./(Z1 + Z2); %parallel R + L with C
Pin = 400;
V = sqrt(abs(Z)*Pin)*sqrt(2);

figure
plot(freq/1E6,V), hold on, grid on, xlabel('Freq (MHz)'),ylabel('Peak Voltage of Vc (V)'),title(['Peak Votlage vs. Frequency - Power In = ',num2str(Pin),'W'])
legend('B = 1T','B = 1.25T','B = 1.5T')

%% Plot Resonances

for j=1:length(C)
    varied_L = L-(L*.1):.01*L:L+(L*.1);
    xL = L-(L*.1):(L+(L*.1))/length(V):L+(L*.1);
    Z1 = R + i.*res_omega(j).*xL;
    Z2 = 1./(i.*res_omega(j)*C(j));
    Z = (Z1.*Z2)./(Z1 + Z2); %parallel R + L with C
    Pin = 400;
    V_new = sqrt(abs(Z)*Pin)*sqrt(2);
    
    figure
    p1 = plot(xL*1E6,V_new), hold on, 
    grid on, 
    xlabel('Varied Inductance (uH)'),
    ylabel('Peak Voltage of Vc (V)'),
    title(['Inductance Drifting - Power In = ',num2str(Pin),'W', ' B = ',num2str(B(j)), 'T'])
    
    % plot %1 points
    x1 = L-(L*.01);
    x1 = x1*1E6;
    y1=get(gca,'ylim')
    hold on
    p2 = plot([x1 x1],y1, '--r')
    x1 = L+(L*.01);
    x1 = x1*1E6;
    y1=get(gca,'ylim')
    hold on
    p7 = plot([x1 x1],y1, '--r')

    % plot %2 points
    x1 = L-(L*.02);
    x1 = x1*1E6;
    y1=get(gca,'ylim')
    hold on
    p3 = plot([x1 x1],y1, '--c')
    x1 = L+(L*.02);
    x1 = x1*1E6;
    y1=get(gca,'ylim')
    hold on
    p8 = plot([x1 x1],y1, '--c')

    % plot %3 points
    x1 = L-(L*.03);
    x1 = x1*1E6;
    y1=get(gca,'ylim')
    hold on
    p4 = plot([x1 x1],y1, '--m')
    x1 = L+(L*.03);
    x1 = x1*1E6;
    y1=get(gca,'ylim')
    hold on
    p9 = plot([x1 x1],y1, '--m')

    % plot %4 points
    x1 = L-(L*.04);
    x1 = x1*1E6;
    y1=get(gca,'ylim')
    hold on
    p5 = plot([x1 x1],y1, '--g')
    x1 = L+(L*.04);
    x1 = x1*1E6;
    y1=get(gca,'ylim')
    hold on
    p10 = plot([x1 x1],y1, '--g')

    % plot %5 points
    x1 = L-(L*.05);
    x1 = x1*1E6;
    y1=get(gca,'ylim')
    hold on
    p6 = plot([x1 x1],y1, '--b')
    x1 = L+(L*.05);
    x1 = x1*1E6;
    y1=get(gca,'ylim')
    hold on
    p11 = plot([x1 x1],y1, '--b')
    legend([p1 p2 p3 p4 p5 p6],'Inductance Drift Curve', '1% variation', '2% variation', '3% variation', '4% variation', '5% variation')
end

%% Inductance variance over capacitance
%want freq to remain constant
for j=1:length(C)
    varied_L = L-(L*.1):.01*L:L+(L*.1);
    set_omega = res_omega(j);      %resonant frew shouldn't change
    c = ones(1,length(varied_L));
    for i=1:length(varied_L)
        cap = 1/(varied_L(i)*set_omega^2);
        c(i) = cap;
    end

    % xC = C-(C*.1):(C+(C*.1))/length(V):C+(C*.1);
    figure
    p1 = plot(varied_L*1E6, c*1E12)
    xlabel('Varied Inductance (uH)'),
    ylabel('Capacitance (pF)')
    title(['Capacitive Correction vs. Inductive Drift',' - B = ',num2str(B(j)), 'T'])

    x1 = L;
    x1 = x1*1E6;
    y1=get(gca,'ylim')
    hold on
    p2 = plot([x1 x1],y1, '-r')

    x1 = L+(0.01*L);
    x1 = x1*1E6;
    y1=get(gca,'ylim')
    hold on
    p3 = plot([x1 x1],y1, '--g')

    x1 = L-(0.01*L);
    x1 = x1*1E6;
    y1=get(gca,'ylim')
    hold on
    plot([x1 x1],y1, '--g')

    x1 = L+(0.09*L);
    x1 = x1*1E6;
    y1=get(gca,'ylim')
    hold on
    p4 = plot([x1 x1],y1, '--c')

    x1 = L-(0.09*L);
    x1 = x1*1E6;
    y1=get(gca,'ylim')
    hold on
    plot([x1 x1],y1, '--c')
    legend([p1 p2 p3 p4],'Inductance Drift Curve', 'Resonance Target Line', '1% variation in L', '9% variation in L')
end

%% Calculate corresponding variables for resonance - Varied L 
res_freq = f_cyc;
res_omega = res_freq * (2*pi);
C = 71E-12;     %Cap of DEE
R = 1;
L = 1./(res_omega.^2.*C);

freq = [5E6:1E2:15E6]'; %sweep from 0 to 20MHz
omega = freq*2*pi;

Z1 = R + i.*omega.*L;
Z2 = 1./(i.*omega*C);
Z = (Z1.*Z2)./(Z1 + Z2); %parallel R + L with C
Pin = 400;
V = sqrt(abs(Z)*Pin)*sqrt(2);

figure
plot(freq/1E6,V), hold on, grid on, xlabel('Freq (MHz)'),ylabel('Peak Voltage of Vc (V)'),title(['Peak Votlage vs. Frequency - Power In = ',num2str(Pin),'W'])
legend('B = 1T','B = 1.25T','B = 1.5T')