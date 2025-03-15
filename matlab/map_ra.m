%initialisation des paramètres
beta1=0.7;
alpha=1;
nu=1;
beta2=0.3;
beta3=0.03;
beta4=0.05;
gamma1=0.01;
gamma2=0.01;
ksi=1;
ki=0.5;
ks=0.3;
omega=0.01;
N=10000;
A=0.095;
Ic=0.25;

%coefficient de diffusion
ds=0.00001;
di=0.1;
dr=0.1;

%discretisation temporelle
T=10;
dt=0.001;
Nt=T/dt;
t=0:dt:T;



mapra=imread("../images/population.png");
mapra2=rgb2gray(imread("../images/populationra_moins.png"));
mapra=rgb2gray(mapra);
figure(1);
imshow(mapra);
figure(2);
imshow(mapra2);

h=1;
size1=size(mapra2);
J2=size1(1);
J1=size1(2);
J=J1*J2;
x=0:h:J2-h;
y=0:h:J1-h;
Sx=zeros(size(mapra2));
Ix=zeros(size(mapra2));
Rx=zeros(size(mapra2));
S=zeros(1,Nt+1);
I=zeros(1,Nt+1);
R=zeros(1,Nt+1);
newSx=zeros(J,1);
newIx=zeros(J,1);
newRx=zeros(J,1);


% Attribution des valeurs selon les conditions spécifiées
Sx(mapra2 >= 248 | mapra2 <= 31) = 0;
Sx(mapra2 > 31 & mapra2 < 248) = 400000 ./ double(mapra2(mapra2 > 31 & mapra2 < 248));

% Affichage de la matrice S
disp(Sx);
% 
figure(3);
im3d=surf(Sx);
set(im3d,'LineStyle','none');
set(gca, 'XDir', 'reverse');
zlim([5, max(Sx(:))])

% Coordonnées du foyer initial d'infectés
foyer_x = 210:220; % Coordonnées x du foyer initial
foyer_y = 120:130; % Coordonnées y du foyer initial

% Définir le foyer initial d'infectés
Ix(foyer_x,foyer_y) = Sx(foyer_x,foyer_y);
Sx(foyer_x,foyer_y) = 0;
Sx=reshape(Sx,J,1);
Ix=reshape(Ix,J,1);
Rx=reshape(Rx,J,1);
S(1)=sum(Sx);
I(1)=sum(Ix);
R(1)=sum(Rx);

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


%figure(3); clf;
%image(x,y,255*reshape(Rx,J2,J1));

% Schema implicite Crank-Nicolson implicit scheme
As =(speye(J) - dt/h^2*ds/2*L0);
Ai=(speye(J) - dt/h^2*di/2*L0);
Ar=(speye(J) - dt/h^2*dr/2*L0);

%boucles
for tt=1:Nt-1
    drawnow;
    %modele simple
    newSx=As\(Sx+(-beta1*alpha*nu*Sx.*Ix-beta2*Sx.*Rx+ds*1/h^2*L0*Sx)*dt);
    newSx(sideleft)=newSx(sideleft+J2);
    newSx(sideright) = newSx(sideright-J2);
    newSx(sidetop) = newSx(sidetop+1);
    newSx(sidebottom) = newSx(sidebottom-1);

    newIx=Ai\(Ix+(beta1*alpha*nu*Sx.*Ix-gamma2*Ix-beta3*Ix.*Rx+di*1/h^2*L0*Ix)*dt);
    newIx(sideleft) = newIx(sideleft+J2);
    newIx(sideright) = newIx(sideright-J2);
    newIx(sidetop) = newIx(sidetop+1);
    newIx(sidebottom) = newIx(sidebottom-1);

    newRx=Ar\(Rx+(gamma2*Ix+beta2*Sx.*Rx+beta3*Ix.*Rx+dr*1/h^2*L0*Rx)*dt);
    newRx(sideleft) = newRx(sideleft+J2);
    newRx(sideright) = newRx(sideright-J2);
    newRx(sidetop) = newRx(sidetop+1);
    newRx(sidebottom) = newRx(sidebottom-1);
    
    %modele avec effet Allee
    %newS=(-beta1*S(tt)*I(tt)*(1-S(tt)/I(tt))*(S(tt)/A-1)-beta2*S(tt)*R(tt))*dt;
    %newI=(beta1*S(tt)*I(tt)*(1-S(tt)/I(tt))*(S(tt)/A-1)-gamma2*I(tt)-beta3*I(tt)*R(tt))*dt;
    %newR=(gamma2*I(tt)+beta2*S(tt)*R(tt)+beta3*I(tt)*R(tt))*dt;
    
    %newSx=As\(Sx+(-beta1*alpha*nu*Sx.*Ix*(1-1/K*I(tt))*(1/A*I(tt)-1)-beta2*Sx.*Rx+omega*Rx+ds*1/h^2*L0*Sx)*dt);
    %newSx=As\(Sx+(-beta1*alpha*nu*Sx.*Ix.*(1-1/K*Ix).*(1/A*Ix-1)-beta2*Sx.*Rx+omega*Rx+ds*1/h^2*L0*Sx)*dt);

    %modele SIRS 
    %newS=(-beta1*S(tt)*I(tt)-beta2*S(tt)*R(tt)+beta4*R(tt))*dt;
    %newI=(beta1*S(tt)*I(tt)-gamma2*I(tt)-beta3*I(tt)*R(tt))*dt;
    %newR=(gamma2*I(tt)+beta2*S(tt)*R(tt)+beta3*I(tt)*R(tt)-beta4*R(tt))*dt;
    % newSx=As\(Sx+(-beta1*Sx.*Ix-beta2*Sx.*Rx+beta4*Rx+ds*1/h^2*L0*Sx)*dt);
    % newSx(sideleft)=newSx(sideleft+J2);
    % newSx(sideright) = newSx(sideright-J2);
    % newSx(sidetop) = newSx(sidetop+1);
    % newSx(sidebottom) = newSx(sidebottom-1);
    % 
    % newIx=Ai\(Ix+(beta1*Sx.*Ix-gamma2*Ix-beta3*Ix.*Rx+di*1/h^2*L0*Ix)*dt);
    % newIx(sideleft) = newIx(sideleft+J2);
    % newIx(sideright) = newIx(sideright-J2);
    % newIx(sidetop) = newIx(sidetop+1);
    % newIx(sidebottom) = newIx(sidebottom-1);
    % 
    % newRx=Ar\(Rx+(gamma2*Ix+beta2*Sx.*Rx+beta3*Ix.*Rx-beta4*Rx+dr*1/h^2*L0*Rx)*dt);
    % newRx(sideleft) = newRx(sideleft+J2);
    % newRx(sideright) = newRx(sideright-J2);
    % newRx(sidetop) = newRx(sidetop+1);
    % newRx(sidebottom) = newRx(sidebottom-1);

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
    
    %avec crédibilité
    % K=0.9;
    % if I(tt)<Ic
    %     newSx=As\(Sx+(-beta1*alpha*nu*Sx.*Ix.*(1-1/K*Ix).*(1/A*Ix-1)-beta2*Sx.*Rx+omega*Rx+ds*1/h^2*L0*Sx)*dt);
    %     newSx(sideleft)=newSx(sideleft+J2);
    %     newSx(sideright) = newSx(sideright-J2);
    %     newSx(sidetop) = newSx(sidetop+1);
    %     newSx(sidebottom) = newSx(sidebottom-1);
    % 
    %     newIx=Ai\(Ix+(beta1*alpha*nu*Sx.*Ix.*(1-1/K*Ix).*(1/A*Ix-1)-gamma2*Ix-beta3*Ix.*Rx+di*1/h^2*L0*Ix)*dt);
    %     newIx(sideleft) = newIx(sideleft+J2);
    %     newIx(sideright) = newIx(sideright-J2);
    %     newIx(sidetop) = newIx(sidetop+1);
    %     newIx(sidebottom) = newIx(sidebottom-1);
    % 
    %     newRx=Ar\(Rx+(gamma2*Ix+beta2*Sx.*Rx+beta3*Ix.*Rx-omega*Rx+dr*1/h^2*L0*Rx)*dt);
    %     newRx(sideleft) = newRx(sideleft+J2);
    %     newRx(sideright) = newRx(sideright-J2);
    %     newRx(sidetop) = newRx(sidetop+1);
    %     newRx(sidebottom) = newRx(sidebottom-1);
    % else
    %     newSx=As\(Sx+(-beta1*alpha*nu*Sx.*Ix.*(1-1/K*Ix).*(1/A*Ix-1)-beta2*Sx.*Rx+omega*Rx+ds*1/h^2*L0*Sx)*dt);
    %     newSx(sideleft)=newSx(sideleft+J2);
    %     newSx(sideright) = newSx(sideright-J2);
    %     newSx(sidetop) = newSx(sidetop+1);
    %     newSx(sidebottom) = newSx(sidebottom-1);
    % 
    %     newIx=Ai\(Ix+(beta1*alpha*nu*Sx.*Ix.*(1-1/K*Ix).*(1/A*Ix-1)-gamma2*Ix-beta3*Ix.*Rx-k*Ix./(1+ksi*Ix)+di*1/h^2*L0*Ix)*dt);
    %     newIx(sideleft) = newIx(sideleft+J2);
    %     newIx(sideright) = newIx(sideright-J2);
    %     newIx(sidetop) = newIx(sidetop+1);
    %     newIx(sidebottom) = newIx(sidebottom-1);
    % 
    %     newRx=Ar\(Rx+(gamma2*Ix+beta2*Sx.*Rx+beta3*Ix.*Rx+k*Ix./(1+ksi*Ix)-omega*Rx+dr*1/h^2*L0*Rx)*dt);
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
    figure(6);
    image(x,y,255*reshape(Ix,J2,J1));
    %figure(3);
    %image(x,y,255*reshape(Rx,J2,J1));

    figure(4);
    im3d=surf(reshape(Ix,J2,J1));
    set(im3d,'LineStyle','none');

    figure(5);
    im3d=surf(reshape(Sx,J2,J1));
    set(im3d,'LineStyle','none');
    zlim([5, max(Sx(:))])
    %if S(tt+1)<0
    %    break
    %end
    %if I(tt+1)<0
    %    break
    %end
end