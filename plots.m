% ACFD reading and plotting data from c++ code
Y(1)=0;Y(102)=1;X(1)=0;X(102)=1;
Y(2:101)=0.005:0.01:0.995;
X(2:101)=0.005:0.01:0.995;
figure(1)
plot(uvdata_upwind(1,:),Y,'k','Linewidth',1.5);
hold on
plot(uvdata_quick(1,:),Y,'r','Linewidth',1.5);
plot(URE100CHOI(:,1),URE100CHOI(:,2),'o','Markeredgecolor','blue');
plot(UGhiaetalRE100(:,1),UGhiaetalRE100(:,2),'s','Markeredgecolor','g','Markersize',8);
ylabel('Y/L','Fontweight' ,'bold');
xlabel('U/Uo','Fontweight' ,'bold');
title({'X velocity variation along center vertical line' '(Lid driven cavity Re=100)'});
legend('upwind','quick','Choi etal','Ghia etal');
ylim([0 1.01]);
hold off

figure(2)
plot(X,uvdata_upwind(2,:),'k','Linewidth',1.5);
hold on
plot(X,uvdata_quick(2,:),'r','Linewidth',1.5);
plot(VchoietalRE100(:,1),VchoietalRE100(:,2),'o','Markeredgecolor','blue');
plot(VRE100ghia(:,1),VRE100ghia(:,2),'s','Markeredgecolor','g','Markersize',8);
xlabel('X/L','Fontweight' ,'bold');
ylabel('V/Uo','Fontweight' ,'bold');
title({'Y velocity variation along center horizontal line' '(Lid driven cavity Re=100)'});
legend('upwind','%quick','Choi etal','Ghia etal');
xlim([-0.01 1.01]);
hold off

figure(3)
contourf(X,Y,udata_upwind,30);
xlabel('X/L','Fontweight' ,'bold','Fontsize',12);
ylabel('Y/L','Fontweight' ,'bold','Fontsize',12);
title({'\fontsize{12}U velocity -- UPWIND' '(lid driven cavity Re=100)'});
colorbar;

figure(4)
contourf(X,Y,udata_quick,30);
xlabel('X/L','Fontweight' ,'bold','Fontsize',12);
ylabel('Y/L','Fontweight' ,'bold','Fontsize',12);
title({'\fontsize{12}U velocity -- QUICK' '(lid driven cavity Re=100)'});
colorbar;


figure(5)
contourf(X,Y,vdata_upwind,30);
xlabel('X/L','Fontweight' ,'bold','Fontsize',12);
ylabel('Y/L','Fontweight' ,'bold','Fontsize',12);
title({'\fontsize{12}V velocity -- UPWIND' '(lid driven cavity Re=100)'});
colorbar;

figure(6)
contourf(X,Y,vdata_quick,30);
xlabel('X/L','Fontweight' ,'bold','Fontsize',12);
ylabel('Y/L','Fontweight' ,'bold','Fontsize',12);
title({'\fontsize{12}V velocity -- QUICK' '(lid driven cavity Re=100)'});
colorbar;

uu = (udata_upwind.^2 + vdata_upwind.^2).^0.5;
vv = (udata_quick.^2 + vdata_quick.^2).^0.5;


figure(7)
contourf(X,Y,uu,30);
xlabel('X/L','Fontweight' ,'bold','Fontsize',12);
ylabel('Y/L','Fontweight' ,'bold','Fontsize',12);
title({'\fontsize{12}total velocity -- UPWIND' '(lid driven cavity Re=100)'});
colorbar;

figure(8)
contourf(X,Y,vv,30);
xlabel('X/L','Fontweight' ,'bold','Fontsize',12);
ylabel('Y/L','Fontweight' ,'bold','Fontsize',12);
title({'\fontsize{12}total velocity -- Quick' '(lid driven cavity Re=100)'});
colorbar;