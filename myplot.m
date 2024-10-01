
figure('position',[100 100 1400 800]);

tiledlayout(3,4,'tilespacing','compact','padding','compact');

nexttile; hold on ; grid on
plot(sigmadata(:,1), sigmadata(:,2)*r2d,'b','LineWidth',2)
text(time(end)*0.7, max(sigmadata(:,2)*r2d*0.7), ...
    ['$t_f = \ $', num2str(time(end),4), '$\ s$'],...
    'fontsize',14,'fontname','times','Interpreter','latex')
xlabel('$t \ (s)$','Interpreter','latex')
ylabel('$\sigma \ (deg)$','Interpreter','latex')
set(gca,'fontsize',14,'fontname','times')

nexttile; hold on ; grid on
plot(V, h,'b','LineWidth',2)
yline(0,':r','LineWidth',2)
text(V(end)*1.3, h(1)*0.6, {['$h(t_f) = \ $', num2str(round(h(end),2)), '$\ km$'] ...
    ['$V(t_f) = \ $', num2str(round(V(end)*1e3,0)), '$\ m/s$']},...
    'fontsize',14,'fontname','times','Interpreter','latex')
xlabel('$V \ (km/s)$','Interpreter','latex')
ylabel('$h \ (km)$','Interpreter','latex')
set(gca,'fontsize',14,'fontname','times')

nexttile; hold on ; grid on
plot(DR, h,'b','LineWidth',2)
yline(0,':r','LineWidth',2)
xlabel('$Downrange \ (km)$','Interpreter','latex')
ylabel('$h \ (km)$','Interpreter','latex')
set(gca,'fontsize',14,'fontname','times')

nexttile; hold on ; grid on
plot(theta, phi,'b','LineWidth',2)
scatter(thetatgt*r2d, phitgt*r2d,100, 'r*')
text(mean(theta), mean(phi), ...
    {['$R_{go}(t_f) = \ $'] [num2str(Rgof, 3), '$\ km$']},...
    'fontsize',14,'fontname','times','Interpreter','latex')
xlabel('$\theta \ (deg)$','Interpreter','latex')
ylabel('$\phi \ (deg)$','Interpreter','latex')
set(gca,'fontsize',14,'fontname','times')

nexttile; hold on ; grid on
plot(time, h,'b','LineWidth',2)
yline(0,':r','LineWidth',2)
xlabel('$t \ (s)$','Interpreter','latex')
ylabel('$h \ (km)$','Interpreter','latex')
set(gca,'fontsize',14,'fontname','times')

nexttile; hold on ; grid on
plot(time, V,'b','LineWidth',2)
xlabel('$t \ (s)$','Interpreter','latex')
ylabel('$V \ (km/s)$','Interpreter','latex')
set(gca,'fontsize',14,'fontname','times')

nexttile; hold on ; grid on
plot(time, gamma,'b','LineWidth',2)
xlabel('$t \ (s)$','Interpreter','latex')
ylabel('$\gamma \ (deg)$','Interpreter','latex')
set(gca,'fontsize',14,'fontname','times')

nexttile; hold on ; grid on
plot(time, psi,'b','LineWidth',2)
xlabel('$t \ (s)$','Interpreter','latex')
ylabel('$\psi \ (deg)$','Interpreter','latex')
set(gca,'fontsize',14,'fontname','times')

nexttile; hold on ; grid on
plot(time, CR,'b','LineWidth',2)
xlabel('$t \ (s)$','Interpreter','latex')
ylabel('$Crossrange \ (km)$','Interpreter','latex')
set(gca,'fontsize',14,'fontname','times')

nexttile; hold on ; grid on
plot(time, A,'b','LineWidth',2)
yline(Amax,':r','LineWidth',2)
ylim([0 Amax*1.1])
xlabel('$t \ (s)$','Interpreter','latex')
ylabel('$A \ (g)$','Interpreter','latex')
set(gca,'fontsize',14,'fontname','times')

nexttile; hold on ; grid on
plot(time, q*1e-3,'b','LineWidth',2)
yline(qmax,':r','LineWidth',2)
ylim([0 qmax*1.1])
xlabel('$t \ (s)$','Interpreter','latex')
ylabel('$q \ (kPa)$','Interpreter','latex')
set(gca,'fontsize',14,'fontname','times')

nexttile; hold on ; grid on
plot(time, qdot,'b','LineWidth',2)
yline(Qdotmax,':r','LineWidth',2)
ylim([0 Qdotmax*1.1])
xlabel('$t \ (s)$','Interpreter','latex')
ylabel('$\dot{Q} \ (kW/m^2)$','Interpreter','latex')
set(gca,'fontsize',14,'fontname','times')