beta1=0.8;
alpha=1;
nu=1;
beta2=0.9;
beta3=0.2;
beta4=0.05;
gamma1=0.2;
gamma2=0.2;

[x1, x2, x3] = meshgrid(0:1:5, 0:1:5,0:1:5);
x1dot = -beta1*alpha*nu*x1.*x2-beta2*x1.*x3; %Note the use of .* and .^
x3dot = gamma2*x2+beta2*x1.*x3+beta3*x2.*x3;
x2dot =beta1*alpha*nu*x1.*x2-gamma2*x2-beta3*x2.*x3;
quiver3(x1,x2,x3,x1dot, x2dot,x3dot);
xlabel('S');ylabel('I');zlabel('R');


% odefun = @(t, x) [-beta1*alpha*nu*x(1).*x(2)-beta2*x(1).*x(3);
%     beta1*alpha*nu*x(1).*x(2)-gamma2*x(2)-beta3*x(2).*x(3);gamma2*x(2)+beta2*x(1).*x(3)+beta3*x(2).*x(3)];
% 
% ICs=[1000, 1, 0]; 
% t=0:.001:50;
% 
% % dD1/dt=-D1/tau;
% % dD2/dt=- D1/tau-Alpha*D2-Beta*h;
% % dh/dt=gamma*D2-Delta*h;
% 
% ICs=[100, 1, 0]; 
% SYS = @(t, X)([-X(1)/tau; 
%               -X(1)/tau-Alpha*X(2)-Beta*X(3); 
%               -Gamma*X(2)-Delta*X(3)]);
% [time, fOUT]=ode45(odefun, t, ICs, [ ]);
% figure(2);
% plot3(fOUT(:,1),fOUT(:,2),fOUT(:,3))
% xlabel('S(t)'), ylabel('I(t)'), zlabel('R(t)');
% title('Solution'), axis tight;
% for x0=0:10:100
%     for y0=0:10:100
%         for z0=0:10:100
%             ICs=[x0, y0, z0]; 
%             SYS = @(t, X)([-X(1)/tau; 
%               -X(1)/tau-Alpha*X(2)-Beta*X(3); 
%               -Gamma*X(2)-Delta*X(3)]);
%             [time, fOUT]=ode45(odefun, t, ICs, [ ]);
%             figure(2); hold on ;
%             plot3(fOUT(:,1),fOUT(:,2),fOUT(:,3))
%         end
%     end
% end