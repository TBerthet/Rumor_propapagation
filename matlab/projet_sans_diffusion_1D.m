%initialisation des paramètres
beta1=0.09;
beta2=0.001;
beta3=0.03;
beta4=0.02;
alpha=0.001;
gamma1=0.001;
gamma2=0.001;
ksi=0.9;
N=10;
A=N/3;

%coefficient de diffusion
ds=0.1;
di=0.1;
dr=0.1;

%discretisation temporelle
T=20;
Nt=1000;
dt=T/Nt;
t=0:dt:T;

%definition de l'espace en 1d
h=0.1;
x0=0;
x1=100;
x=x0:h:x1;
J=length(x);

%variables dynamiques
S=zeros(1,Nt+1);
I=zeros(1,Nt+1);
R=zeros(1,Nt+1);
newS=zeros(1,Nt+1);
newI=zeros(1,Nt+1);
newR=zeros(1,Nt+1);

%conditions initiales
S(1)=990;
I(1)=10;
R(1)=0;
%S(1)=(beta1(2*beta1-2*beta2+beta3)*N+gamma1*(beta1-beta2))/beta1*(beta1-beta2+beta3);
%I(1)=(beta1*beta2*N+gamma*beta2)/beta1*(beta1-beta2+beta3);
%condition de Neumann et laplacien
 L = sparse(1:J,1:J,-2); % matrice creuse, compacte en memoire/sparse matrix, compact in memory, with -2 on the diagonal
 L = spdiags(ones(J,2),[-1 1],L); % forme la matrice tridiagonale/fill in the off-diagonals 
 L(1,:) = 0; % v(x0) est donne par les conditions au bord/v(x0) will be set by the boundary conditions
 L(J,:) = 0; % w(x1) est donne par les conditions au bord/w(x1) will be set by the boundary conditions

%boucles
for tt=1:Nt-1
    %modele simple
    newS=(-beta1*S(tt)*I(tt)-beta2*S(tt)*R(tt))*dt;
    newI=(beta1*S(tt)*I(tt)-gamma2*I(tt)-beta3*I(tt)*R(tt))*dt;
    newR=(gamma2*I(tt)+beta2*S(tt)*R(tt)+beta3*I(tt)*R(tt))*dt;
    
    %modele avec effet Allee
    %newS=(-beta1*S(tt)*I(tt)*(1-S(tt)/I(tt))*(S(tt)/A-1)-beta2*S(tt)*R(tt))*dt;
    %newI=(beta1*S(tt)*I(tt)*(1-S(tt)/I(tt))*(S(tt)/A-1)-gamma2*I(tt)-beta3*I(tt)*R(tt))*dt;
    %newR=(gamma2*I(tt)+beta2*S(tt)*R(tt)+beta3*I(tt)*R(tt))*dt;

    %modele SIRS 
    %newS=(-beta1*S(tt)*I(tt)-beta2*S(tt)*R(tt)+beta4*R(tt))*dt;
    %newI=(beta1*S(tt)*I(tt)-gamma2*I(tt)-beta3*I(tt)*R(tt))*dt;
    %newR=(gamma2*I(tt)+beta2*S(tt)*R(tt)+beta3*I(tt)*R(tt)-beta4*R(tt))*dt;

    %modele avec incidence dependant de S^2
    %newS=(-beta1*S(tt)^2*I(tt)-beta2*S(tt)*R(tt))*dt;
    %newI=(beta1*S(tt)^2*I(tt)-gamma2*I(tt)-beta3*I(tt)*R(tt))*dt;
    %newR=(gamma2*I(tt)+beta2*S(tt)*R(tt)+beta3*I(tt)*R(tt))*dt;

    %modele avec taux d'incidence non monotone
    %peut faire de jolies choses
    %newS=(-beta1*S(tt)*I(tt)/(1+alpha*I(tt)^2)-beta2*S(tt)*R(tt))*dt;
    %newI=(beta1*S(tt)*I(tt)/(1+alpha*I(tt)^2)-gamma2*I(tt)-beta3*I(tt)*R(tt))*dt;
    %newR=(gamma2*I(tt)+beta2*S(tt)*R(tt)+beta3*I(tt)*R(tt))*dt;

    S(tt+1)=S(tt)+newS;
    I(tt+1)=I(tt)+newI;
    R(tt+1)=R(tt)+newR;
    %if S(tt+1)<0
    %    break
    %end
    %if I(tt+1)<0
    %    break
    %end
end

    
%plot
plot(t,S,t,I,t,R);
grid on;
xlabel("Time t");
ylabel("Number of individuals");
legend('S','I','R');

f=@(t,Y)[-beta1*Y(1)*Y(2)-beta2*Y(1)*(N-Y(1)-Y(2))+beta4*(N-Y(1)-Y(2));beta1*Y(1)*Y(2)-gamma2*Y(2)^2-beta3*Y(2)*(N-Y(1)-Y(2))];
y1=linspace(0,N,20);
y2=linspace(0,N,20);
[x,y]=meshgrid(y1,y2);
S=zeros(size(x));
I=zeros(size(x));
t=0;
for i=1:numel(x)
    Yprime=f(t,[x(i); y(i)]);
    S(i)=Yprime(1);
    I(i)=Yprime(2);
end
quiver(x,y,S,I,'r');figure(gcf)
xlabel('S');
ylabel('I');
axis tight equal;
hold on
for y20 = 0:5
    for x20=5
        [ts,ys] = ode45(f,[0,400],[x20;y20]);
        plot(ys(:,1),ys(:,2))
        plot(ys(1,1),ys(1,2),'bo') % starting point
        plot(ys(end,1),ys(end,2),'ks') % ending point
    end
end
hold off
