%initialisation des paramètres
beta1=0.5;
alpha=1;
nu=1;
beta2=0.3;
beta3=0.03;
beta4=0.05;
gamma1=0.01;
gamma2=0.01;
ksi=1;
k=0.5;
omega=0.5;
A=6;
K=40;
Ic=0.25;
%coefficient de diffusion
ds=0.1;
di=0.1;
dr=0.1;

%discretisation temporelle
T=30;
dt=0.01;
Nt=T/dt;
t=0:dt:T;

%definition de l'espace en 1d
h=0.1;
x0=0;
x1=10;
y0=0;
y1=10;
x=x0:h:(x1-h);
y=y0:h:(y1-h);
[X,Y]=meshgrid(x,y);
J1=length(x);
J2=length(y);
J=J1*J2;

%variables dynamiques
S=zeros(1,Nt+1);
I=zeros(1,Nt+1);
R=zeros(1,Nt+1);
P=zeros(J,1);
Sx=zeros(J,1);
Ix=zeros(J,1);
Rx=zeros(J,1);
newSx=zeros(J,1);
newIx=zeros(J,1);

newRx=zeros(J,1);

%conditions initiales
% Initialisation des conditions initiales pour Sx et Ix
S_initial = 20; % Proportion initiale de personnes susceptibles
I_initial = 30; % Proportion initiale de personnes infectées

% Coordonnées du foyer initial d'infectés
foyer_x = 50; % Coordonnées x du foyer initial
foyer_y = 50; % Coordonnées y du foyer initial

% Initialisation de Sx et Ix
Sx = ones(J, 1)*50; % Initialisation de Sx avec la valeur initiale

% Définir le foyer initial d'infectés
Ix(foyer_x + foyer_y*J1) = I_initial;
Sx(foyer_x + foyer_y*J1) = S_initial;
% Sx=0.7+ 0.1*(-1 + 2*rand(J,1));
% Ix=.3 + 0.1*(-1 + 2*rand(J,1));
Rx(1)=1;
Sx(1)=0;
N=sum(Sx)+sum(Ix)+sum(Rx);
S(1)=sum(Sx)/N;
I(1)=sum(Ix)/N;
R(1)=sum(Rx)/N;

% discretisation du Laplacien avec conditions aux bord de Neumann
% no flux, i.e. du/dx = 0 en x0 et en x1
% Neumann, no flux boundary conditions

% Le Laplacien discretise est une matrice L de taille JxJ
% Discretized Laplacian 

cornertopleft = 1;
cornerbottomleft = J2;
cornertopright = (J1-1)*J2+1;
cornerbottomright = J;
sideleft = 2:J2-1;
sidetop = J2+1:J2:J2*(J1-2)+1;
sidebottom = 2*J2:J2:J2*(J1-1);
sideright = J2*(J1-1)+2:J-1;
side = [cornertopleft, cornertopright, cornerbottomleft, cornerbottomright, ...
    sideleft, sidetop, sidebottom, sideright];
interieur = setdiff(1:J, side);

% interieur
L0 = sparse(interieur,interieur,-4,J,J); % matrice creuse, compacte en memoire
L0 = L0 + sparse(interieur,interieur+1,1,J,J);
L0 = L0 + sparse(interieur,interieur-1,1,J,J);
L0 = L0 + sparse(interieur,interieur+J2,1,J,J);
L0 = L0 + sparse(interieur,interieur-J2,1,J,J);

figure(1); clf;
image(x,y,255*reshape(Ix,J2,J1));

%figure(3); clf;
%image(x,y,255*reshape(Rx,J2,J1));

% Schema implicite Crank-Nicolson implicit scheme
As = (speye(J) - dt/h^2*ds/2*L0);
Ai=(speye(J) - dt/h^2*di/2*L0);
Ar=(speye(J) - dt/h^2*dr/2*L0);
%boucles

temps=0;
savet=0:0.1:T;
for tt=1:Nt
    drawnow;
    %modele simple
    % newSx=As\(Sx+(-beta1*alpha*nu*Sx.*Ix-beta2*Sx.*Rx+ds*1/h^2*L0*Sx)*dt);
    % newSx(sideleft)=newSx(sideleft+J2);
    % newSx(sideright) = newSx(sideright-J2);
    % newSx(sidetop) = newSx(sidetop+1);
    % newSx(sidebottom) = newSx(sidebottom-1);
    % 
    % newIx=Ai\(Ix+(beta1*alpha*nu*Sx.*Ix-gamma2*Ix-beta3*Ix.*Rx+di*1/h^2*L0*Ix)*dt);
    % newIx(sideleft) = newIx(sideleft+J2);
    % newIx(sideright) = newIx(sideright-J2);
    % newIx(sidetop) = newIx(sidetop+1);
    % newIx(sidebottom) = newIx(sidebottom-1);
    % 
    % newRx=Ar\(Rx+(gamma2*Ix+beta2*Sx.*Rx+beta3*Ix.*Rx+dr*1/h^2*L0*Rx)*dt);
    % newRx(sideleft) = newRx(sideleft+J2);
    % newRx(sideright) = newRx(sideright-J2);
    % newRx(sidetop) = newRx(sidetop+1);
    % newRx(sidebottom) = newRx(sidebottom-1);
    
    % modele avec effet Allee
    % newSx=As\(Sx+(-beta1*alpha*nu*Sx.*Ix.*(1-Ix./K).*(Ix./A-1)-beta2*Sx.*Rx+ds*1/h^2*L0*Sx)*dt);
    % newSx(sideleft)=newSx(sideleft+J2);
    % newSx(sideright) = newSx(sideright-J2);
    % newSx(sidetop) = newSx(sidetop+1);
    % newSx(sidebottom) = newSx(sidebottom-1);
    % 
    % newIx=Ai\(Ix+(beta1*alpha*nu*Sx.*Ix.*(1-Ix./K).*(Ix./A-1)-gamma2*Ix-beta3*Ix.*Rx+di*1/h^2*L0*Ix)*dt);
    % newIx(sideleft) = newIx(sideleft+J2);
    % newIx(sideright) = newIx(sideright-J2);
    % newIx(sidetop) = newIx(sidetop+1);
    % newIx(sidebottom) = newIx(sidebottom-1);
    % 
    % newRx=Ar\(Rx+(gamma2*Ix+beta2*Sx.*Rx+beta3*Ix.*Rx+dr*1/h^2*L0*Rx)*dt);
    % newRx(sideleft) = newRx(sideleft+J2);
    % newRx(sideright) = newRx(sideright-J2);
    % newRx(sidetop) = newRx(sidetop+1);
    % newRx(sidebottom) = newRx(sidebottom-1);

    %modele SIRS 
    newSx=As\(Sx+(-beta1*alpha*nu*Sx.*Ix-beta2*Sx.*Rx+omega*Rx+ds*1/h^2*L0*Sx)*dt);
    newSx(sideleft)=newSx(sideleft+J2);
    newSx(sideright) = newSx(sideright-J2);
    newSx(sidetop) = newSx(sidetop+1);
    newSx(sidebottom) = newSx(sidebottom-1);

    newIx=Ai\(Ix+(beta1*alpha*nu*Sx.*Ix-gamma2*Ix-beta3*Ix.*Rx+di*1/h^2*L0*Ix)*dt);
    newIx(sideleft) = newIx(sideleft+J2);
    newIx(sideright) = newIx(sideright-J2);
    newIx(sidetop) = newIx(sidetop+1);
    newIx(sidebottom) = newIx(sidebottom-1);

    newRx=Ar\(Rx+(gamma2*Ix+beta2*Sx.*Rx+beta3*Ix.*Rx-omega*Rx+dr*1/h^2*L0*Rx)*dt);
    newRx(sideleft) = newRx(sideleft+J2);
    newRx(sideright) = newRx(sideright-J2);
    newRx(sidetop) = newRx(sidetop+1);
    newRx(sidebottom) = newRx(sidebottom-1);

    %modele avec incidence dependant de S^2
    %newS=(-beta1*S(tt)^2*I(tt)-beta2*S(tt)*R(tt))*dt;
    %newI=(beta1*S(tt)^2*I(tt)-gamma2*I(tt)-beta3*I(tt)*R(tt))*dt;
    %newR=(gamma2*I(tt)+beta2*S(tt)*R(tt)+beta3*I(tt)*R(tt))*dt;

    %modele avec taux d'incidence non monotone
    %newSx=Sx+(-beta1*Sx.*Ix./(1+alpha*Ix.^2)-beta2*Sx.*Rx+ds*1/h^2*L*Sx)*dt;
    %newIx=Ix+(beta1*Sx.*Ix./(1+alpha*Ix.^2)-gamma2*Ix-beta3*Ix.*Rx+di*1/h^2*L*Ix)*dt;
    %newRx=Rx+(gamma2*Ix+beta2*Sx.*Rx+beta3*Ix.*Rx+dr*1/h^2*L*Rx)*dt;
    % newSx=As\(Sx+(-beta1*Sx.*Ix./(1+alpha*Ix.^2)-beta2*Sx.*Rx+ds*1/h^2*L0*Sx)*dt);
    % newSx(sideleft)=newSx(sideleft+J2);
    % newSx(sideright) = newSx(sideright-J2);
    % newSx(sidetop) = newSx(sidetop+1);
    % newSx(sidebottom) = newSx(sidebottom-1);
    % 
    % newIx=Ai\(Ix+(beta1*Sx.*Ix./(1+alpha*Ix.^2)-gamma2*Ix-beta3*Ix.*Rx+di*1/h^2*L0*Ix)*dt);
    % newIx(sideleft) = newIx(sideleft+J2);
    % newIx(sideright) = newIx(sideright-J2);
    % newIx(sidetop) = newIx(sidetop+1);
    % newIx(sidebottom) = newIx(sidebottom-1);
    % 
    % newRx=Ar\(Rx+(gamma2*Ix+beta2*Sx.*Rx+beta3*Ix.*Rx+dr*1/h^2*L0*Rx)*dt);
    % newRx(sideleft) = newRx(sideleft+J2);
    % newRx(sideright) = newRx(sideright-J2);
    % newRx(sidetop) = newRx(sidetop+1);
    % newRx(sidebottom) = newRx(sidebottom-1);
    
    % avec contrôle
    % K=0.9;
    % if I(tt)<Ic
    %     newSx=As\(Sx+(-beta1*alpha*nu*Sx.*Ix-beta2*Sx.*Rx+ds*1/h^2*L0*Sx)*dt);
    %     newSx(sideleft)=newSx(sideleft+J2);
    %     newSx(sideright) = newSx(sideright-J2);
    %     newSx(sidetop) = newSx(sidetop+1);
    %     newSx(sidebottom) = newSx(sidebottom-1);
    % 
    %     newIx=Ai\(Ix+(beta1*alpha*nu*Sx.*Ix-gamma2*Ix-beta3*Ix.*Rx+di*1/h^2*L0*Ix)*dt);
    %     newIx(sideleft) = newIx(sideleft+J2);
    %     newIx(sideright) = newIx(sideright-J2);
    %     newIx(sidetop) = newIx(sidetop+1);
    %     newIx(sidebottom) = newIx(sidebottom-1);
    % 
    %     newRx=Ar\(Rx+(gamma2*Ix+beta2*Sx.*Rx+beta3*Ix.*Rx+dr*1/h^2*L0*Rx)*dt);
    %     newRx(sideleft) = newRx(sideleft+J2);
    %     newRx(sideright) = newRx(sideright-J2);
    %     newRx(sidetop) = newRx(sidetop+1);
    %     newRx(sidebottom) = newRx(sidebottom-1);
    % else
    %     newSx=As\(Sx+(-beta1*alpha*nu*Sx.*Ix-beta2*Sx.*Rx-k*Sx./(1+ksi*Sx)+ds*1/h^2*L0*Sx)*dt);
    %     newSx(sideleft)=newSx(sideleft+J2);
    %     newSx(sideright) = newSx(sideright-J2);
    %     newSx(sidetop) = newSx(sidetop+1);
    %     newSx(sidebottom) = newSx(sidebottom-1);
    % 
    %     newIx=Ai\(Ix+(beta1*alpha*nu*Sx.*Ix-gamma2*Ix-beta3*Ix.*Rx-k*Ix./(1+ksi*Ix)+di*1/h^2*L0*Ix)*dt);
    %     newIx(sideleft) = newIx(sideleft+J2);
    %     newIx(sideright) = newIx(sideright-J2);
    %     newIx(sidetop) = newIx(sidetop+1);
    %     newIx(sidebottom) = newIx(sidebottom-1);
    % 
    %     newRx=Ar\(Rx+(gamma2*Ix+beta2*Sx.*Rx+beta3*Ix.*Rx+k*Ix./(1+ksi*Ix)+k*Sx./(1+ksi*Sx)+dr*1/h^2*L0*Rx)*dt);
    %     newRx(sideleft) = newRx(sideleft+J2);
    %     newRx(sideright) = newRx(sideright-J2);
    %     newRx(sidetop) = newRx(sidetop+1);
    %     newRx(sidebottom) = newRx(sidebottom-1);
    % end
    newSx(1)=newSx(2);
    newIx(1)=newIx(2);
    newRx(1)=newRx(2);
    Sx=newSx;
    Ix=newIx;
    Rx=newRx;
    figure(1);
    image(x,y,255*reshape(Ix./(Ix+Sx+Rx),J2,J1));
    title(['t = ' num2str(temps)]);
    xlabel("x");ylabel("y");
    % nom_fichier = ['plot_2d_', num2str(temps), '.png'];
    % chemin_complet = fullfile('C:\Users\berth\OneDrive\Bureau\4bim\S1\biomaths4\projetbm4\rumor_propagation\picture\oublie\o_0.001', nom_fichier);
    % saveas(gcf,chemin_complet);

    % Fermer la figure actuelle pour éviter l'accumulation
    %figure(3);
    %image(x,y,255*reshape(Rx,J2,J1));
    % figure(4);
    % surf(reshape(Ix./(Ix+Sx+Rx),J2,J1));
    % zlabel("Densité de I");
    % title(['t = ' num2str(temps)]);
    % xlabel("x");
    % ylabel("y");
    S(tt+1)=sum(Sx)/N;
    I(tt+1)=sum(Ix)/N;
    R(tt+1)=sum(Rx)/N;
    % P=(1-1/K*Ix).*(1/A*Ix-1);
    % figure(5);
    % plot(Ix,P);
    %if S(tt+1)<0
    %    break
    %end
    %if I(tt+1)<0
    %    break
    %end
    temps=temps+dt;
end

    
%plot

figure(2);
plot(t,S,t,I,t,R);
grid on;
xlabel("Time t");
ylabel("Density of individuals");
legend('S','I','R');