function [] = projective_geometry_demo(sz)
% A "pedagogic" demo of projective geometry showing some mathematical
% objects projections. Tune the projection parameters (view point,
% projection matrix) to test your own and change the output.
%
% Author & support nicolas.douillet (at) free.fr, 2007-2020.


[X,Y] = meshgrid(0:sz-1);
Z = zeros(1,4);

% Part I : classic 2D display (top view)

% Chessboard 2D display
figure;
set(gcf,'Color',[1 1 1]);
subplot(121);

for i = 1:sz   
    for j = 1:sz      
        if ~mod(i+j,2)            
            clr = [0 0 0];            
        elseif mod(i+j,2)            
            clr = [1 1 0];            
        end        
        if i < sz && j < sz            
            fill3([X(i,j) X(i,j+1) X(i+1,j+1) X(i+1,j)],...
                  [Y(i,j) Y(i,j+1) Y(i+1,j+1) Y(i+1,j)],Z,clr), hold on;            
        end        
    end    
end

axis tight, axis square;
grid on;
view(2);
title('Chessboard, circle, parabole, and hyperbole, standard geometry, top view','FontSize',16);


% Denominators
D = X+Y+1;

% Draw hyperbole
hyper_nb_pts = 2e3;
H = compute_hyperbole(hyper_nb_pts,sz);
Xh = H(1,:);
Yh = H(2,:);
plot(Xh,Yh,'Color',[0 0 1],'Linewidth',2), hold on;

% Draw circle / ellipse
Rho = 1;
a = 1;
b = 1;
circle_nb_pts = 1e3;
angle_step = 2*pi/circle_nb_pts;
theta = 0:angle_step:2*pi-angle_step;
Xe = Rho*(1 + a*cos(theta));
Ye = Rho*(1 + b*sin(theta));
plot(Xe,Ye,'Color',[1 0 0],'Linewidth',2), hold on;

% Draw parabola
para_nb_pts = 7e2;
P = compute_parabola(para_nb_pts);
plot(P(1,1:para_nb_pts),P(2,1:para_nb_pts),'Color',[0 1 1],'LineWidth',2), hold on;
plot(P(1,para_nb_pts+1:end),P(2,para_nb_pts+1:end),'Color',[0 1 1],'LineWidth',2), hold on;
axis([0 sz-1 0 sz-1]);


% Part II : projections and projective geometry items display (from [2;0] point of view)

% Compute lines equation in projective geometry
% 
% X’ = (4X+2)/(X+Y+1)
% Y’ = (X+Y)/(X+Y+1) 
% D = X+Y+1 % denominator

M = [4,0;1,1]; % TRANSFORMATION MATRIX
C = [2;0]; % VIEW POINT

X = X(:)';
Y = Y(:)';
Cxy = cat(1,X,Y);
D = D(:)';

% Chessboard projection
Cxyp = (M*Cxy+repmat(C,[1,sz^2]))./repmat(D,[2,1]);

Xp = reshape(Cxyp(1,:)',[sz,sz]);
Yp = reshape(Cxyp(2,:)',[sz,sz]);

subplot(122);

for i = 1:sz   
    for j = 1:sz      
        if ~mod(i+j,2)            
            clr = [0 0 0];            
        elseif mod(i+j,2)            
            clr = [1 1 0];            
        end        
        if i < sz && j < sz            
            fill3([Xp(i,j) Xp(i,j+1) Xp(i+1,j+1) Xp(i+1,j)],...
                  [Yp(i,j) Yp(i,j+1) Yp(i+1,j+1) Yp(i+1,j)],Z,clr), hold on;            
        end        
    end    
end

axis tight, axis square;
box on, grid on;
view(2);


% Circle / ellipse projection
Exy = cat(1,Xe,Ye);
Exyp = (M*Exy+repmat(C,[1,circle_nb_pts]))./repmat(Xe+Ye+1,[2,1]);
plot(Exyp(1,:),Exyp(2,:),'Color',[1 0 0],'Linewidth',2), hold on;


% Hyperbole projection
Hxy = cat(1,Xh,Yh);
Hpxy = (M*Hxy+repmat(C,[1,hyper_nb_pts]))./repmat(Xh+Yh+1,[2,1]);
plot(Hpxy(1,:),Hpxy(2,:),'Color',[0 0 1],'Linewidth',2), hold on;

% Hyperbole negative coordinates
Hnxy = (M*(-Hxy)+repmat(C,[1,hyper_nb_pts]))./repmat(-Xh-Yh+1,[2,1]);
plot(Hnxy(1,:),Hnxy(2,:),'Color',[0 0 1],'Linewidth',2), hold on;


% Parabole projection
Pxyp = (M*P+repmat(C,[1,2*para_nb_pts]))./repmat(1+sum(P,1),[2,1]);
plot(Pxyp(1,1:para_nb_pts),Pxyp(2,1:para_nb_pts),'Color',[0 1 1],'Linewidth',2), hold on;
plot(Pxyp(1,para_nb_pts+1:end),Pxyp(2,para_nb_pts+1:end),'Color',[0 1 1],'Linewidth',2), hold on;

title('Chessboard, circle, parabole, and hyperbole, projective geometry view','FontSize',16);

end


function [P] = compute_parabola(para_nb_pts)
% Parabole iterative construction technique 
% 
% Author & support nicolas.douillet (at) free.fr, 2007-2020.


k = 0.5*sqrt(2);
R = k*[1, -1; 1, 1]; % pi/4 2D rotation matrix

% Parabola summit : (0.5,0.5)

t = 1:para_nb_pts;
P = zeros(2,size(t,2));

P(1,:) = ((0.1*(t-1)).^2)/2/k + k; % before rotation : F(1,0)
P(2,:) = 0.1*(t-1);
P = R*P;
P = cat(2,P,flipud(P));


end


function [H] = compute_hyperbole(hyper_nb_pts,sz)
% Hyperbole iterative construction technique 
% 
% Author & support nicolas.douillet (at) free.fr, 2007-2020.

I = eye(2);
J = [0, 1;-1, 0];  % pi/2 2D rotation matrix;
C = [0,0.5;0.5,0];
theta_h = 2e-2;

E = (cosh(theta_h))*I+(sinh(theta_h))*J*C; % or E = expm(theta_h*J*C); % matrix exponential
H = zeros(2,hyper_nb_pts);
H(:,1) = [1/sz;sz];  % hyperbole inital point

for k = 2:hyper_nb_pts
    H(:,k) = E*H(:,k-1);
end


end