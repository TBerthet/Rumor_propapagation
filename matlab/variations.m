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
ds=0.1;
di=0.1;
dr=0.1;

%discretisation temporelle
T=300;
dt=0.4;
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
Sx=zeros(J,1);
Ix=zeros(J,1);
Rx=zeros(J,1);
newSx=zeros(J,1);
newIx=zeros(J,1);
newRx=zeros(J,1);


%conditions initiales
Sx(x>2)=1;
Ix(x<=2)=1;
S(1)=sum(Sx)/J;
I(1)=sum(Ix)/J;
R(1)=sum(Rx)/J;

Sx1=Sx;
Sx2=Sx;
Sx3=Sx;
Sx4=Sx;
Sx5=Sx;
Ix1=Ix;
Ix2=Ix;
Ix3=Ix;
Ix4=Ix;
Ix5=Ix;
Rx1=Rx;
Rx2=Rx;
Rx3=Rx;
Rx4=Rx;
Rx5=Rx;

S1=S;
S2=S;
S3=S;
S4=S;
S5=S;
I1=I;
I2=I;
I3=I;
I4=I;
I5=I;
R1=R;
R2=R;
R3=R;
R4=R;
R5=R;

Stx1=Sx1;
Itx1=Ix1;
Rtx1=Rx1;

Stx3=Sx3;
Itx3=Ix3;
Rtx3=Rx3;

Stx5=Sx5;
Itx5=Ix5;
Rtx5=Rx5;
%condition de Neumann et laplacien
 L = sparse(1:J,1:J,-2); % matrice creuse, compacte en memoire/sparse matrix, compact in memory, with -2 on the diagonal
 L = spdiags(ones(J,2),[-1 1],L); % forme la matrice tridiagonale/fill in the off-diagonals 
 L(1,:) = 0; % v(x0) est donne par les conditions au bord/v(x0) will be set by the boundary conditions
 L(J,:) = 0; % w(x1) est donne par les conditions au bord/w(x1) will be set by the boundary conditions
%boucles
beta31=-0.5;
beta32=-0.05;
beta33=-0.01;
beta34=-0.005;
beta35=0;
for tt=1:Nt
    %modele simple
    newSx1=Sx1+(-beta1*alpha*nu*Sx1.*Ix1-beta2*Sx1.*Rx1+ds*1/h^2*L*Sx1)*dt;
    newIx1=Ix1+(beta1*alpha*nu*Sx1.*Ix1-gamma2*Ix1-beta31*Ix1.*Rx1+di*1/h^2*L*Ix1)*dt;
    newRx1=Rx1+(gamma2*Ix1+beta2*Sx1.*Rx1+beta31*Ix1.*Rx1+dr*1/h^2*L*Rx1)*dt;
   
    newSx1(1)=newSx1(2);
    newIx1(1)=newIx1(2);
    newRx1(1)=newRx1(2);

    Sx1=newSx1;
    Ix1=newIx1;
    Rx1=newRx1;
    Stx1=cat(2,Stx1,Sx1);
    Itx1=cat(2,Itx1,Ix1);
    Rtx1=cat(2,Rtx1,Rx1);
    figure(1);
    plot(x,Sx1,x,Ix1,x,Rx1,x,Sx1+Ix1+Rx1);
    grid on;
    xlabel("x");
    ylabel("Proportion of individuals");
    legend('S','I','R','somme');
    S1(tt+1)=sum(Sx1)/J;
    I1(tt+1)=sum(Ix1)/J;
    R1(tt+1)=sum(Rx1)/J;
    
    newSx2=Sx2+(-beta1*alpha*nu*Sx2.*Ix2-beta2*Sx2.*Rx2+ds*1/h^2*L*Sx2)*dt;
    newIx2=Ix2+(beta1*alpha*nu*Sx2.*Ix2-gamma2*Ix2-beta32*Ix2.*Rx2+di*1/h^2*L*Ix2)*dt;
    newRx2=Rx2+(gamma2*Ix2+beta2*Sx2.*Rx2+beta32*Ix2.*Rx2+dr*1/h^2*L*Rx2)*dt;
   
    newSx2(1)=newSx2(2);
    newIx2(1)=newIx2(2);
    newRx2(1)=newRx2(2);

    Sx2=newSx2;
    Ix2=newIx2;
    Rx2=newRx2;
    S2(tt+1)=sum(Sx2)/J;
    I2(tt+1)=sum(Ix2)/J;
    R2(tt+1)=sum(Rx2)/J;
    
    %I3
    newSx3=Sx3+(-beta1*alpha*nu*Sx3.*Ix3-beta2*Sx3.*Rx3+ds*1/h^2*L*Sx3)*dt;
    newIx3=Ix3+(beta1*alpha*nu*Sx3.*Ix3-gamma2*Ix3-beta33*Ix3.*Rx3+di*1/h^2*L*Ix3)*dt;
    newRx3=Rx3+(gamma2*Ix3+beta2*Sx3.*Rx3+beta33*Ix3.*Rx3+dr*1/h^2*L*Rx3)*dt;
   
    newSx3(1)=newSx3(2);
    newIx3(1)=newIx3(2);
    newRx3(1)=newRx3(2);

    Sx3=newSx3;
    Ix3=newIx3;
    Rx3=newRx3;
    Stx3=cat(2,Stx3,Sx3);
    Itx3=cat(2,Itx3,Ix3);
    Rtx3=cat(2,Rtx3,Rx3);
    S3(tt+1)=sum(Sx3)/J;
    I3(tt+1)=sum(Ix3)/J;
    R3(tt+1)=sum(Rx3)/J;

    %I4
    newSx4=Sx4+(-beta1*alpha*nu*Sx4.*Ix4-beta2*Sx4.*Rx4+ds*1/h^2*L*Sx4)*dt;
    newIx4=Ix4+(beta1*alpha*nu*Sx4.*Ix4-gamma2*Ix4-beta34*Ix4.*Rx4+di*1/h^2*L*Ix4)*dt;
    newRx4=Rx4+(gamma2*Ix4+beta2*Sx4.*Rx4+beta34*Ix4.*Rx4+dr*1/h^2*L*Rx4)*dt;
   
    newSx4(1)=newSx4(2);
    newIx4(1)=newIx4(2);
    newRx4(1)=newRx4(2);

    Sx4=newSx4;
    Ix4=newIx4;
    Rx4=newRx4;
    S4(tt+1)=sum(Sx4)/J;
    I4(tt+1)=sum(Ix4)/J;
    R4(tt+1)=sum(Rx4)/J;

    %I5
    newSx5=Sx5+(-beta1*alpha*nu*Sx5.*Ix5-beta2*Sx5.*Rx5+ds*1/h^2*L*Sx5)*dt;
    newIx5=Ix5+(beta1*alpha*nu*Sx5.*Ix5-gamma2*Ix5-beta35*Ix5.*Rx5+di*1/h^2*L*Ix5)*dt;
    newRx5=Rx5+(gamma2*Ix5+beta2*Sx5.*Rx5+beta35*Ix5.*Rx5+dr*1/h^2*L*Rx5)*dt;
   
    newSx5(1)=newSx5(2);
    newIx5(1)=newIx5(2);
    newRx5(1)=newRx5(2);

    Sx5=newSx5;
    Ix5=newIx5;
    Rx5=newRx5;
    Stx5=cat(2,Stx5,Sx5);
    Itx5=cat(2,Itx5,Ix5);
    Rtx5=cat(2,Rtx5,Rx5);
    S5(tt+1)=sum(Sx5)/J;
    I5(tt+1)=sum(Ix5)/J;
    R5(tt+1)=sum(Rx5)/J;

end

    
%plot
figure(3);
plot(t,I1,t,I2,t,I3,t,I4,t,I5);
grid on;
xlabel("Time t");
ylabel("density of individuals I");
legend('\beta_3=-0.5','\beta_3=-0.05','\beta_3=-0.01', '\beta_3=-0.005', '\beta_3=0')
figure(2);
h1=surf(Itx1)
set(h1,'LineStyle','none');
xlabel("time t");
ylabel("x");
zlabel("I")
figure(4);
h3=surf(Itx3)
set(h3,'LineStyle','none');
xlabel("time t");
ylabel("x");
zlabel("I")
figure(5);
h5=surf(Itx5)
set(h5,'LineStyle','none');
xlabel("time t");
ylabel("x");
zlabel("I")