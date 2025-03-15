%initialisation des paramètres
beta1=0.5;
alpha=1;
nu=1;
beta2=0.3;
beta3=0.03;
beta4=0.05;
gamma1=0.2;
gamma2=0.002;
gamma3=0.01;
ksi=0.9;
N=100;
A=N/2;

%coefficient de diffusion
ds=1;
di=0.1;
dr=0.1;

%discretisation temporelle
T=30;
dt=0.05;
Nt=T/dt;
t=0:dt:T;

%definition de l'espace en 1d
h=1;
x0=0;
x1=50;
x=x0:h:x1;
J=length(x);

%variables dynamiques
S=zeros(1,Nt+1);
I=zeros(1,Nt+1);
R=zeros(1,Nt+1);
Sx=50*ones(J,1);
Ix=zeros(J,1);
Rx=zeros(J,1);
newSx=zeros(J,1);
newIx=zeros(J,1);
newRx=zeros(J,1);


%conditions initiales
Ix(x<=2)=5;
Sx(x<=2)=45;
S(1)=sum(Sx)/sum(Sx+Rx+Ix);
I(1)=sum(Ix)/sum(Sx+Rx+Ix);
R(1)=sum(Rx)/sum(Sx+Rx+Ix);
Stx=Sx;
Itx=Ix./(Sx+Ix+Rx);
Rtx=Rx;
%condition de Neumann et laplacien
 L = sparse(1:J,1:J,-2); % matrice creuse, compacte en memoire/sparse matrix, compact in memory, with -2 on the diagonal
 L = spdiags(ones(J,2),[-1 1],L); % forme la matrice tridiagonale/fill in the off-diagonals 
 L(1,:) = 0; % v(x0) est donne par les conditions au bord/v(x0) will be set by the boundary conditions
 L(J,:) = 0; % w(x1) est donne par les conditions au bord/w(x1) will be set by the boundary conditions
%boucles

for tt=1:Nt
    %modele simple
    newSx=Sx+(-beta1*alpha*nu*Sx.*Ix-beta2*Sx.*Rx+ds*1/h^2*L*Sx)*dt;
    newIx=Ix+(beta1*alpha*nu*Sx.*Ix-gamma2*Ix-beta3*Ix.*Rx+di*1/h^2*L*Ix)*dt;
    newRx=Rx+(gamma2*Ix+beta2*Sx.*Rx+beta3*Ix.*Rx+dr*1/h^2*L*Rx)*dt;
    
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
    % newSx=Sx+(-beta1*Sx.*Ix./(1+alpha*Ix.^2)-beta2*Sx.*Rx+ds*1/h^2*L*Sx)*dt;
    % newIx=Ix+(beta1*Sx.*Ix./(1+alpha*Ix.^2)-gamma2*Ix-beta3*Ix.*Rx+gamma3*Ix.*Rx+di*1/h^2*L*Ix)*dt;
    % newRx=Rx+(gamma2*Ix+beta2*Sx.*Rx+beta3*Ix.*Rx-gamma3*Ix.*Rx+dr*1/h^2*L*Rx)*dt;
    newSx(1)=newSx(2);
    newIx(1)=newIx(2);
    newRx(1)=newRx(2);

    Sx=newSx;
    Ix=newIx;
    Rx=newRx;
    Stx=cat(2,Stx,Sx);
    Itx=cat(2,Itx,Ix./(Sx+Ix+Rx));
    Rtx=cat(2,Rtx,Rx);
    figure(1);
    plot(x,Sx,x,Ix,x,Rx,x,Sx+Ix+Rx);
    grid on;
    xlabel("x");
    ylabel("Proportion of individuals");
    legend('S','I','R','somme');
    %figure(2);
    %plot(x,Sx);
    S(tt+1)=sum(Sx)/sum(Sx+Rx+Ix);
    I(tt+1)=sum(Ix)/sum(Sx+Rx+Ix);
    R(tt+1)=sum(Rx)/sum(Sx+Rx+Ix);
    
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
ylabel("proportions of individuals");
legend('S','I','R');
figure(3);
plot(t,I);
grid on;
xlabel("Time t");
ylabel("density of individuals I");
figure(2);

h=surf(Itx)
set(h,'LineStyle','none');
xlabel("time t");
ylabel("x");
zlabel("I")