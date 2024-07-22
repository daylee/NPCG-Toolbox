% myplot.m

figure('position',[100 60 1300 800]);

tiledlayout(3,4,'tilespacing','compact','padding','compact');

nexttile; hold on ; grid on
plot(sigmadata(:,1), sigmadata(:,2)*r2d,'b','linewidth',2)
text(time(end)*0.7,  siglmt*r2d*0.9, ...
    ['t_f = ', num2str(time(end),4), ' s'],'fontsize',17,'fontname','times')
xlabel('time (s)')
ylabel('\sigma (deg)')
%ylim([-siglmt siglmt]*r2d)
set(gca,'fontsize',16,'fontname','times')

nexttile; hold on ; grid on
plot(V, h,'b','linewidth',2)
yline(0,':r','LineWidth',3)
text(V(end)*1.2, h(1)*0.6, {['h(t_f) = ', num2str(round(h(end),2)), ' km'] ...
['V(t_f) = ', num2str(round(V(end)*1e3,2)), ' m/s']},'fontsize',17,'fontname','times')
xlabel('V (km/s)')
ylabel('h (km)')
set(gca,'fontsize',16,'fontname','times')

nexttile; hold on ; grid on
plot(DR, h,'b','linewidth',2)
yline(0,':r','LineWidth',3)
xlabel('Downrange (km)')
ylabel('h (km)')
set(gca,'fontsize',16,'fontname','times')

nexttile; hold on ; grid on
plot(theta, phi,'b','linewidth',2)
scatter(thetatgt*r2d, phitgt*r2d,100, 'r*')
xlabel('\theta (deg)')
ylabel('\phi (deg)')
set(gca,'fontsize',16,'fontname','times')

nexttile; hold on ; grid on
plot(time, h,'b','linewidth',2)
yline(0,':r','LineWidth',3)
xlabel('time (s)')
ylabel('h (km)')
set(gca,'fontsize',16,'fontname','times')

nexttile; hold on ; grid on
plot(time, V,'b','linewidth',2)
xlabel('time (s)')
ylabel('V (km/s)')
set(gca,'fontsize',16,'fontname','times')

nexttile; hold on ; grid on
plot(time, gamma,'b','linewidth',2)
xlabel('time (s)')
ylabel('\gamma (deg)')
set(gca,'fontsize',16,'fontname','times')
 
nexttile; hold on ; grid on
plot(time, psi,'b','linewidth',2)
xlabel('time (s)')
ylabel('\psi (deg)')
set(gca,'fontsize',16,'fontname','times')

nexttile; hold on ; grid on
plot(time, CR,'b','linewidth',2)
xlabel('time (s)')
ylabel('Crossrange (km)')
yline(0,':r','LineWidth',3)
set(gca,'fontsize',16,'fontname','times')

nexttile; hold on ; grid on
plot(time, A,'b','linewidth',2)
yline(Amax,':r','LineWidth',3)
ylim([0 Amax*1.1])
xlabel('time (s)')
ylabel('A (g)')
set(gca,'fontsize',16,'fontname','times')

nexttile; hold on ; grid on
plot(time, q*1e-3,'b','linewidth',2)
yline(qmax,':r','LineWidth',3)
ylim([0 qmax*1.1])
xlabel('time (s)')
ylabel('q (kPa)')
set(gca,'fontsize',16,'fontname','times')

nexttile; hold on ; grid on
plot(time, qdot,'b','linewidth',2)
yline(Qdotmax,':r','LineWidth',3)
ylim([0 Qdotmax*1.1])
xlabel('time (s)')
ylabel('Qdot (kW/m^2)')
set(gca,'fontsize',16,'fontname','times')

