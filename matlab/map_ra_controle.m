%initialisation des paramètres
beta1=0.5;
alpha=1;
nu=0.6;
beta2=0.2;
beta3=0.005;
beta4=0.05;
gamma1=0.01;
gamma2=1/30;

omega=0.01;
A=0.095;
Ic=0.25;
ksi=1;
k=0.5;

%coefficient de diffusion
ds=0.1;
di=0.1;
dr=0.1;

%discretisation temporelle
T=10;
dt=0.0005;
Nt=T/dt;
t=0:dt:T;



mapra=imread("../images/population.png");
mapra2=rgb2gray(imread("../images/populationra_moins.png"));
mapra=rgb2gray(mapra);
figure(1);
imshow(mapra);
figure(2);
imshow(mapra2);

h=10;
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
% figure(3);
% im3d=surf(Sx);
% set(im3d,'LineStyle','none');
% zlim([5, max(Sx(:))])

% Coordonnées du foyer initial d'infectés
foyer_x = 165:167; % Coordonnées x du foyer initial
foyer_y = 157:159; % Coordonnées y du foyer initial

% Définir le foyer initial d'infectés
Ix(foyer_x,foyer_y) = 40;
Rx(foyer_x,foyer_y) = 10;
Sx(foyer_x,foyer_y) = Sx(foyer_x,foyer_y)-50;
Sx=reshape(Sx,J,1);
Ix=reshape(Ix,J,1);
Rx=reshape(Rx,J,1);
N=sum(Sx)+sum(Rx)+sum(Ix);
S(1)=sum(Sx)/N;
I(1)=sum(Ix)/N;
R(1)=sum(Rx)/N;

% % Condition periodiques
% % periodic conditions
% LX = sparse(1:J,1:J,-2); % matrice creuse, compacte en memoire/sparse matrix, compact in memory, with -2 on the diagonal
% LY = sparse(1:J,1:J,-2); % matrice creuse, compacte en memoire/sparse matrix, compact in memory, with -2 on the diagonal
% upperleftcorner = 1;
% lowerleftcorner = J2;
% upperrightcorner = J2*(J1-1)+1;
% lowerrightcorner = J;
% leftborder = 2:J2-1;
% upperborder = J2+1:J2:J2*(J1-2)+1;
% lowerborder = 2*J2:J2:J2*(J1-1);
% rightborder = J2*(J1-1)+2:J-1;
% border = [upperleftcorner, upperrightcorner, lowerleftcorner, lowerrightcorner, ...
%     leftborder, upperborder, lowerborder, rightborder];
% interiorY = setdiff(1:J, [upperborder,lowerborder, ...
%   upperleftcorner,upperrightcorner,lowerleftcorner,lowerrightcorner]);
% 
% [ix,iy] = ind2sub([J2,J1],setdiff(1:J,[leftborder,rightborder, ...
%   upperleftcorner,upperrightcorner,lowerleftcorner,lowerrightcorner]));
% interiorX = sub2ind([J1,J2],iy,ix);
% 
% %% interior
% % On construit deux matrices du Laplacien discretise: LX et LY. LX est le Laplacien
% % discretise en prenant les indices en X d'abord:
% % We build two matrices for the discretized Laplacian: LX and LY. LX
% % is the Laplacian taking indices along the x-axis first.
% % 
% %   1    2 3 4 5 ...   J1
% %   J1+1 ...         2*J1
% %   ...
% %   J1*(J2-1)+1 ... J1*J2
% %
% % LY est le Laplacien discretise en prenant les indices le long de Y d'abord:
% % LY is the disctretized Laplacian taking indices along the y-axis first 
% %  
% %   1 J2+1 ...(J1-1)*J2+1
% %   2 J2+2            ...
% %   3 
% %   ...               ...
% %   J2 2*J2 3*J2 ... J1*J2
% %
% % LX et LY sont des matrices tridiagonales
% % LX and LY are tri-diagonal matrices 
% 
% % interior
% LX = LX + sparse(interiorX,interiorX+1,1,J,J);
% LX = LX + sparse(interiorX,interiorX-1,1,J,J);
% 
% LY = LY + sparse(interiorY,interiorY+1,1,J,J);
% LY = LY + sparse(interiorY,interiorY-1,1,J,J);
% 
% % borders
% LY = LY + sparse(upperborder,upperborder+1,1,J,J);
% LY = LY + sparse(upperborder,upperborder+J2-1,1,J,J);
% 
% [ix,iy] = ind2sub([J2,J1],leftborder);
% leftborderX = sub2ind([J1,J2],iy,ix);
% LX = LX + sparse(leftborderX,leftborderX+1,1,J,J);
% LX = LX + sparse(leftborderX,leftborderX+J1-1,1,J,J);
% 
% LY = LY + sparse(lowerborder,lowerborder-(J2-1),1,J,J);
% LY = LY + sparse(lowerborder,lowerborder-1,1,J,J);
% 
% [ix,iy] = ind2sub([J2,J1],rightborder);
% rightborderX = sub2ind([J1,J2],iy,ix);
% LX = LX + sparse(rightborderX,rightborderX+1,1,J,J);
% LX = LX + sparse(rightborderX,rightborderX+J1-1,1,J,J);
% 
% % corners
% LY(upperleftcorner,upperleftcorner+1) = 1;
% LY(upperleftcorner,upperleftcorner+J2-1) = 1;
% 
% [ix,iy] = ind2sub([J2,J1],upperleftcorner);
% upperleftcornerX = sub2ind([J1,J2],iy,ix);
% LX(upperleftcornerX,upperleftcornerX+1) = 1;
% LX(upperleftcornerX,upperleftcornerX+J1-1) = 1;
% 
% 
% LY(lowerleftcorner,lowerleftcorner-(J2-1)) = 1;
% LY(lowerleftcorner,lowerleftcorner-1) = 1;
% 
% [ix,iy] = ind2sub([J2,J1],lowerleftcorner);
% lowerleftcornerX = sub2ind([J1,J2],iy,ix);
% LX(lowerleftcornerX,lowerleftcornerX+1) = 1;
% LX(lowerleftcornerX,lowerleftcornerX+J1-1) = 1;
% 
% LY(upperrightcorner,upperrightcorner+1) = 1;
% LY(upperrightcorner,upperrightcorner+J2-1) = 1;
% 
% [ix,iy] = ind2sub([J2,J1],upperrightcorner);
% upperrightcornerX = sub2ind([J1,J2],iy,ix);
% LX(upperrightcornerX,upperrightcornerX-(J1-1)) = 1;
% LX(upperrightcornerX,upperrightcornerX-1) = 1;
% 
% LY(lowerrightcorner,lowerrightcorner-(J2-1)) = 1;
% LY(lowerrightcorner,lowerrightcorner-1) = 1;
% 
% [ix,iy] = ind2sub([J2,J1],lowerrightcorner);
% lowerrightcornerX = sub2ind([J1,J2],iy,ix);
% LX(lowerrightcornerX,lowerrightcornerX-(J1-1)) = 1;
% LX(lowerrightcornerX,lowerrightcornerX-1) = 1;

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
interiorY = setdiff(1:J, side);

[ix,iy] = ind2sub([J2,J1],interiorY);
interiorX = sub2ind([J1,J2],iy,ix);

%% test indices Y
% MY = zeros(J2,J1);
% MY(interiorY) = 1;
% spy(MY)
%% test indices X
% MX = zeros(J1,J2);
% MX(interiorX) = 1;
% spy(MX)

%% interior
% On construit deux matrices du Laplacien discretise: LX et LY. LX est le Laplacien
% discretise en prenant les indices en X d'abord:
% We build two matrices for the discretized Laplacian: LX and LY. LX
% is the Laplacian taking indices along the x-axis first.
% 
%   1    2 3 4 5 ...   J1
%   J1+1 ...         2*J1
%   ...
%   J1*(J2-1)+1 ... J1*J2
%
% LY est le Laplacien discretise en prenant les indices le long de Y d'abord:
% LY is the disctretized Laplacian taking indices along the y-axis first 
%  
%   1 J2+1 ...(J1-1)*J2+1
%   2 J2+2            ...
%   3 
%   ...               ...
%   J2 2*J2 3*J ... J1*J2
%
% LX et LY sont des matrices tridiagonales
% LX and LY are tri-diagonal matrices 

LX = sparse(interiorX,interiorX,-2,J,J); % matrice creuse, compacte en memoire
LX = LX + sparse(interiorX,interiorX+1,1,J,J);
LX = LX + sparse(interiorX,interiorX-1,1,J,J);

LY = sparse(interiorY,interiorY,-2,J,J); % matrice creuse, compacte en memoire
LY = LY + sparse(interiorY,interiorY+1,1,J,J);
LY = LY + sparse(interiorY,interiorY-1,1,J,J);

% Schema implicite Crank-Nicolson implicit scheme
% ADI: alternate direction implicit method
AXs = (speye(J) - dt/h^2*ds/2*LX);
AYs = (speye(J) - dt/h^2*ds/2*LY);
AXi = (speye(J) - dt/h^2*di/2*LX);
AYi = (speye(J) - dt/h^2*di/2*LY);
AXr = (speye(J) - dt/h^2*dr/2*LX);
AYr = (speye(J) - dt/h^2*dr/2*LY);

figure(5);
im3d=surf(reshape(Sx,J2,J1));
set(im3d,'LineStyle','none');
zlim([5, max(Sx(:))])

%boucles
temps=0;
for tt=1:Nt-1
    drawnow;
    % Methode ADI:
    % On resout en 2 etapes: d'abord on fait par rapport a X:
    % u12 = u + k/2*r*u*(1-u) + k/2/h^2*D*(LX*u12 + LY*u)
    % ensuite par rapport a Y:
    % u = u12 + k/2*r*u12*(1-u12) + k/2/h^2*D*(LX*u12 + LY*u)
    % L'evaluation se fait a k/2 a chaque etape

    % ADI method: 
    % Two step solver: first we apply diffusion only along the x-axis
    % u12 = u + k/2*r*u*(1-u) + k/2/h^2*D*(LX*u12 + LY*u)
    % then we make a change of index and apply diffusion along the y-axis
    % u = u12 + k/2*r*u12*(1-u12) + k/2/h^2*D*(LX*u12 + LY*u)
    % Each step advances by k/2 time step 
    %modele simple
    if I(tt)<Ic
        bb = (Sx+(-beta1*alpha*nu*Sx.*Ix-beta2*Sx.*Rx+ds*1/h^2*LY*Sx)*dt/2); % En X: on calcule les donnee du systeme implicite AX*u12 = b
                                                       % Along X: evaluate data b for the implicite system AX*u12 = b            
        bb = reshape(reshape(bb,J2,J1)',J,1);            % On convertit les donnees en indices X-dominant
                                                       % Convert data to X-dominant index             
        s12 = AXs\bb;                                    % On resout
                                                       % solve
        LXs12 = LX*s12;                                % On stocke LX*u12
                                                       % Store LX*u12 for later
        s12 = reshape(reshape(s12,J1,J2)',J,1);        % On reconvertit u12 en Y-dominant
                                                       % Convert back u12 to Y-dominant index
        LXs12 = reshape(reshape(LXs12,J1,J2)',J,1);    % on reconvertit LXu12 en Y-dominant
                                                       % Convert LXu12 to Y-dominant
        bb = (Ix+(beta1*alpha*nu*Sx.*Ix-gamma2*Ix-beta3*Ix.*Rx+di*1/h^2*LY*Ix)*dt/2); % En X: on calcule les donnee du systeme implicite AX*u12 = b
                                                       % Along X: evaluate data b for the implicite system AX*u12 = b            
        bb = reshape(reshape(bb,J2,J1)',J,1);            % On convertit les donnees en indices X-dominant
                                                       % Convert data to X-dominant index             
        i12 = AXi\bb;                                    % On resout
                                                       % solve
        LXi12 = LX*i12;                                % On stocke LX*u12
                                                       % Store LX*u12 for later
        i12 = reshape(reshape(i12,J1,J2)',J,1);        % On reconvertit u12 en Y-dominant
                                                       % Convert back u12 to Y-dominant index
        LXi12 = reshape(reshape(LXi12,J1,J2)',J,1);    % on reconvertit LXu12 en Y-dominant
                                                       % Convert LXu12 to Y-dominant
    
        bb = (Rx+(gamma2*Ix+beta2*Sx.*Rx+beta3*Ix.*Rx+dr*1/h^2*LY*Rx)*dt/2); % En X: on calcule les donnee du systeme implicite AX*u12 = b
                                                       % Along X: evaluate data b for the implicite system AX*u12 = b            
        bb = reshape(reshape(bb,J2,J1)',J,1);            % On convertit les donnees en indices X-dominant
                                                       % Convert data to X-dominant index             
        r12 = AXr\bb;                                    % On resout
                                                       % solve
        LXr12 = LX*r12;                                % On stocke LX*u12
                                                       % Store LX*u12 for later
        r12 = reshape(reshape(r12,J1,J2)',J,1);        % On reconvertit u12 en Y-dominant
                                                       % Convert back u12 to Y-dominant index
        LXr12 = reshape(reshape(LXr12,J1,J2)',J,1);    % on reconvertit LXu12 en Y-dominant
                                                       % Convert LXu12 to Y-dominant
    
        newSx = AYs\(s12+(-beta1*alpha*nu*s12.*i12-beta2*s12.*r12+ds*1/h^2*LXs12)*dt/2); % En Y: on resout AY*u = b 
                                                                   % Along Y: solve AY*u = b
    
        newIx = AYi\(i12+(beta1*alpha*nu*s12.*i12-gamma2*i12-beta3*i12.*r12+di*1/h^2*LXi12)*dt/2); % En Y: on resout AY*u = b 
                                                                   % Along Y: solve AY*u = b
    
        newRx=AYr\(r12+(gamma2*i12+beta2*s12.*r12+beta3*i12.*r12+dr*1/h^2*LXr12)*dt/2);
        % end of ADI 
    else
         bb = (Sx+(-beta1*alpha*nu*Sx.*Ix-beta2*Sx.*Rx-k*Sx./(1+ksi*Sx)+ds*1/h^2*LY*Sx)*dt/2); % En X: on calcule les donnee du systeme implicite AX*u12 = b
                                                       % Along X: evaluate data b for the implicite system AX*u12 = b            
        bb = reshape(reshape(bb,J2,J1)',J,1);            % On convertit les donnees en indices X-dominant
                                                       % Convert data to X-dominant index             
        s12 = AXs\bb;                                    % On resout
                                                       % solve
        LXs12 = LX*s12;                                % On stocke LX*u12
                                                       % Store LX*u12 for later
        s12 = reshape(reshape(s12,J1,J2)',J,1);        % On reconvertit u12 en Y-dominant
                                                       % Convert back u12 to Y-dominant index
        LXs12 = reshape(reshape(LXs12,J1,J2)',J,1);    % on reconvertit LXu12 en Y-dominant
                                                       % Convert LXu12 to Y-dominant
        bb = (Ix+(beta1*alpha*nu*Sx.*Ix-gamma2*Ix-beta3*Ix.*Rx-k*Ix./(1+ksi*Ix)+di*1/h^2*LY*Ix)*dt/2); % En X: on calcule les donnee du systeme implicite AX*u12 = b
                                                       % Along X: evaluate data b for the implicite system AX*u12 = b            
        bb = reshape(reshape(bb,J2,J1)',J,1);            % On convertit les donnees en indices X-dominant
                                                       % Convert data to X-dominant index             
        i12 = AXi\bb;                                    % On resout
                                                       % solve
        LXi12 = LX*i12;                                % On stocke LX*u12
                                                       % Store LX*u12 for later
        i12 = reshape(reshape(i12,J1,J2)',J,1);        % On reconvertit u12 en Y-dominant
                                                       % Convert back u12 to Y-dominant index
        LXi12 = reshape(reshape(LXi12,J1,J2)',J,1);    % on reconvertit LXu12 en Y-dominant
                                                       % Convert LXu12 to Y-dominant
    
        bb = (Rx+(gamma2*Ix+beta2*Sx.*Rx+beta3*Ix.*Rx+k*Sx./(1+ksi*Sx)+k*Ix./(1+ksi*Ix)+dr*1/h^2*LY*Rx)*dt/2); % En X: on calcule les donnee du systeme implicite AX*u12 = b
                                                       % Along X: evaluate data b for the implicite system AX*u12 = b            
        bb = reshape(reshape(bb,J2,J1)',J,1);            % On convertit les donnees en indices X-dominant
                                                       % Convert data to X-dominant index             
        r12 = AXr\bb;                                    % On resout
                                                       % solve
        LXr12 = LX*r12;                                % On stocke LX*u12
                                                       % Store LX*u12 for later
        r12 = reshape(reshape(r12,J1,J2)',J,1);        % On reconvertit u12 en Y-dominant
                                                       % Convert back u12 to Y-dominant index
        LXr12 = reshape(reshape(LXr12,J1,J2)',J,1);    % on reconvertit LXu12 en Y-dominant
                                                       % Convert LXu12 to Y-dominant
    
        newSx = AYs\(s12+(-beta1*alpha*nu*s12.*i12-beta2*s12.*r12-k*s12./(1+ksi*s12)+ds*1/h^2*LXs12)*dt/2); % En Y: on resout AY*u = b 
                                                                   % Along Y: solve AY*u = b
    
        newIx = AYi\(i12+(beta1*alpha*nu*s12.*i12-gamma2*i12-beta3*i12.*r12-k*i12./(1+ksi*i12)+di*1/h^2*LXi12)*dt/2); % En Y: on resout AY*u = b 
                                                                   % Along Y: solve AY*u = b
    
        newRx=AYr\(r12+(gamma2*i12+beta2*s12.*r12+beta3*i12.*r12+k*s12./(1+ksi*s12)+k*i12./(1+ksi*i12)+dr*1/h^2*LXr12)*dt/2);
        % end of ADI 
    end
    Sx=newSx;
    Ix=newIx;
    Rx=newRx;
    Sx(sideleft) = Sx(sideleft+J2);              % condition Neumann no flux condition
    Sx(sideright) = Sx(sideright-J2);
    Sx(sidetop) = Sx(sidetop+1);
    Sx(sidebottom) = Sx(sidebottom-1);
    Ix(sideleft) = Ix(sideleft+J2);              % condition Neumann no flux condition
    Ix(sideright) = Ix(sideright-J2);
    Ix(sidetop) = Ix(sidetop+1);
    Ix(sidebottom) = Ix(sidebottom-1);
    Rx(sideleft) = Rx(sideleft+J2);              % condition Neumann no flux condition
    Rx(sideright) = Rx(sideright-J2);
    Rx(sidetop) = Rx(sidetop+1);
    Rx(sidebottom) = Rx(sidebottom-1);
    figure(6);
    image(x,y,255*reshape(Ix,J2,J1));
    title(['t = ' num2str(temps)]);
    xlabel("x");
    ylabel("y");
    % nom_fichier = ['plot_2d_I_', num2str(temps), '.png'];
    % chemin_complet = fullfile('C:\Users\berth\OneDrive\Bureau\4bim\S1\biomaths4\projetbm4\rumor_propagation\picture\contexte\rumeur_nulle', nom_fichier);
    % saveas(gcf,chemin_complet);
    % figure(8);
    % image(x,y,255*reshape(Rx,J2,J1));
    % title(['t = ' num2str(temps)]);
    % xlabel("x");
    % ylabel("y");

    % figure(4);
    % im3d=surf(reshape(Ix,J2,J1));
    % set(im3d,'LineStyle','none');

    % figure(5);
    % im3d=surf(reshape(Sx,J2,J1));
    % set(im3d,'LineStyle','none');
    % zlim([5, max(Sx(:))])
    
    
    % figure(7);
    % image(x,y,255*reshape(Sx+Rx,J2,J1));
    % title(['t = ' num2str(temps)]);
    % xlabel("x");
    % ylabel("y");
    % nom_fichier = ['plot_2d_S_', num2str(temps), '.png'];
    % chemin_complet = fullfile('C:\Users\berth\OneDrive\Bureau\4bim\S1\biomaths4\projetbm4\rumor_propagation\picture\contexte\rumeur_nulle', nom_fichier);
    % saveas(gcf,chemin_complet);
    temps=temps+dt;
    S(tt+1)=sum(Sx)/N;
    I(tt+1)=sum(Ix)/N;
    R(tt+1)=sum(Rx)/N;
    %if S(tt+1)<0
    %    break
    %end
    %if I(tt+1)<0
    %    break
    %end
end