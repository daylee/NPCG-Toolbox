
figure('position',[50 300 1000 650]);
tiledlayout(3,3,'tilespacing','compact','padding','compact');

nexttile; hold on ; grid on
plot(sigmadata(:,1), sigmadata(:,2)*r2d,'b','linewidth',2)
text(time(end)*0.7, min(sigmadata(:,2))*r2d*0.8, ...
    ['t_f = ', num2str(time(end),4), ' s'],'fontsize',17,'fontname','times')
xlabel('t (s)')
ylabel('\sigma (deg)')
% ylim([-siglmt siglmt]*r2d)
set(gca,'fontsize',16,'fontname','times')

nexttile; hold on ; grid on
plot(V, h,'b','linewidth',2)
yline(0,':r','LineWidth',2)
text(V(end)*1.2, h(1)*0.6, {['h(t_f) = ', num2str(round(h(end),2)), ' km'] ...
    ['V(t_f) = ', num2str(round(V(end)*1e3,2)), ' m/s']},'fontsize',17,'fontname','times')
xlabel('V (km/s)')
ylabel('h (km)')
set(gca,'fontsize',16,'fontname','times')

nexttile; hold on ; grid on
plot(time, gamma,'b','linewidth',2)
text(time(end)*0.3, gamma(end)*0.9, ...
    ['\gamma(t_f) = ', num2str(gamma(end),2), ' deg'],'fontsize',17,'fontname','times')
xlabel('t (s)')
ylabel('\gamma (deg)')
set(gca,'fontsize',16,'fontname','times')

nexttile; hold on ; grid on
plot(time, DRtogo,'b','linewidth',2)
text(time(end)*0.5, DRtogo(1)*0.5, ...
    ['R^D_{go}(t_f) = ', num2str(round(abs(DRtogo(end)),2)), ' km'],'fontsize',17,'fontname','times')
xlabel('t (s)')
ylabel('R^D_{go} (km)')
set(gca,'fontsize',16,'fontname','times')

nexttile; hold on ; grid on
plot(time, CR,'b','linewidth',2)
text(time(end)*0.5, min(CR)*0.5, ...
    ['R^C(t_f) = ', num2str(round(abs(CR(end)),2)), ' km'],'fontsize',17,'fontname','times')
xlabel('t (s)')
ylabel('R^C (km)')
set(gca,'fontsize',16,'fontname','times')

nexttile; hold on ; grid on
plot(theta, phi,'b','linewidth',2)
scatter(thetatgt*r2d, phitgt*r2d,100, 'r*')
text( 0.5*(max(theta)+min(theta)), 0.5*(max(phi)+min(phi)),...
    ['R_{go}(t_f) = ', num2str(Rgof, 2), ' km'],'fontsize',17,'fontname','times')
xlabel('\theta (deg)')
ylabel('\phi (deg)')
set(gca,'fontsize',16,'fontname','times')

nexttile; hold on ; grid on
plot(time, A,'b','linewidth',2)
yline(Amax,':r','LineWidth',2)
text(time(end)*0.3, max(A)*1.1, ...
    ['max(A) = ', num2str(max(A),3), ' g'],'fontsize',17,'fontname','times')
ylim([0 Amax*1.1])
xlabel('t (s)')
ylabel('A (g)')
set(gca,'fontsize',16,'fontname','times')

nexttile; hold on ; grid on
plot(time, q*1e-3,'b','linewidth',2)
yline(qmax,':r','LineWidth',2)
text(time(end)*0.3, max(q*1e-3)*1.1, ...
    ['max(q) = ', num2str(max(q*1e-3),3), ' kPa'],'fontsize',17,'fontname','times')
ylim([0 qmax*1.1])
xlabel('t (s)')
ylabel('q (kPa)')
set(gca,'fontsize',16,'fontname','times')

nexttile; hold on ; grid on
plot(time, qdot,'b','linewidth',2)
yline(Qdotmax,':r','LineWidth',2)
text(time(end)*0.2, max(qdot)*1.1, ...
    ['max(Qdot) = ', num2str(max(qdot),4), ' kW/m^2'],'fontsize',17,'fontname','times')
ylim([0 Qdotmax*1.1])
xlabel('t (s)')
ylabel('Qdot (kW/m^2)')
set(gca,'fontsize',16,'fontname','times')

