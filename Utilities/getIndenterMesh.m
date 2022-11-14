function [F,V] = getIndenterMesh(alpha_lim, phi_lim, R, ns, pointSpacing,resampleCurveOpt,interpMethod)
%% Plot indenter flag
plot_indenter = 0;

%% based on HELP_regionTriMesh3D
min_alpha = deg2rad(alpha_lim(1));
max_alpha = deg2rad(alpha_lim(2));
phi = deg2rad(phi_lim);

% Plot settings
fontSize=15;
markerSize1=45;
lineWidth1=4;
faceAlpha=0.5;
% pointSpacing_factor = 10;

%Boundary 1
t = []; r = []; x_temp = []; XY_temp = []; z = []; x = []; y = [];
t=flip(linspace(phi(1),phi(2),round(ns*(phi(2)*R))));
t(end) = [];% t(1) = [];
r=R*ones(size(t));
[z,x_temp] = pol2cart(t,r); z = R-z;
XY_temp = x_temp'*[cos(max_alpha),sin(max_alpha)];
x = XY_temp(:,1); y = XY_temp(:,2);

V1=[x(:) y(:) z(:)];
%Boundary 2
t = []; r = []; x_temp = []; XY_temp = []; z = []; x = []; y = [];

t=linspace(phi(1),phi(2),round(ns*(phi(2)*R)));
t(end) = []; %t(1)  = [];
r=R*ones(size(t));
[z,x_temp] = pol2cart(t,r); z = R-z;
XY_temp = x_temp'*[cos(min_alpha),sin(min_alpha)];
x = XY_temp(:,1); y = XY_temp(:,2);

V2=[x(:) y(:) z(:)];
%Boundary 3
t = []; r = []; x_temp = []; XY_temp = []; z = []; x = []; y = [];

t=linspace(min_alpha,max_alpha,round(ns*(R*(max_alpha-min_alpha))));
t(end) = [];
[x,y] = pol2cart(t,sin(phi(2))*R*ones(size(t)));
z=R*(1-cos(phi(2)))*ones(size(x));
V3=[x(:) y(:) z(:)];
% V = [V1;V2;V3];
V = [V1;V2];
regionCell={V};

% >>>>>>>>>>>>>>>>plot<<<<<<<<<<<<<<<<<<<<<<
if plot_indenter    
    figure;
    hold on;
    scatter3(V1(:,1),V1(:,2),V1(:,3), 'o');
    scatter3(V2(:,1),V2(:,2),V2(:,3), '+');
    scatter3(V3(:,1),V3(:,2),V3(:,3), 'x');
    plot3(V(:,1),V(:,2),V(:,3));
end
% >>>>>>>>>>>>>>>>plot<<<<<<<<<<<<<<<<<<<<<<

%% 
% Meshing the region (See also |regionTriMesh2D|)

%Defining a region and control parameters (See also |regionTriMesh2D|)
% pointSpacing=0.05; %Desired point spacing
% resampleCurveOpt=1; 
% interpMethod='linear'; %or 'natural'

mustPointsBoundary(1,:) = [0 0 0];
mustPointsBoundary(2,:) = [R*cos(min_alpha)*(sin(phi(2))) R*sin(min_alpha)*(sin(phi(2))) R*(1-cos(phi(2)))];
mustPointsBoundary(3,:) = [R*cos(max_alpha)*(sin(phi(2))) R*sin(max_alpha)*(sin(phi(2))) R*(1-cos(phi(2)))];
% mustPointsBoundary = [];
[F,V]=regionTriMesh3D(regionCell,pointSpacing,resampleCurveOpt,interpMethod);


% [~, min_x_ind]=min(V(:,1)); %find the nearest node from the origin (checked with x value)
% V(min_x_ind,:) = [0 0 0]; % relocate the node to the origion.

[V,~] = setCorner(V,[0 0 0]);
[V,I2] = setCorner(V,[R*cos(min_alpha)*(sin(phi(2))) R*sin(min_alpha)*(sin(phi(2))) R*(1-cos(phi(2)))]);
[V,I3] = setCorner(V,[R*cos(max_alpha)*(sin(phi(2))) R*sin(max_alpha)*(sin(phi(2))) R*(1-cos(phi(2)))]);
if I2==I3
    V(end,:) = mustPointsBoundary(2,:);
end


% 
% [~, max_y_ind]=max(V(:,2)); 
% V(max_y_ind,:) = [R*cos(max_alpha) R*sin(max_alpha) R]; % relocate the node to the origion.
% 
% [~, min_y_ind]=min(V(:,2)); 
% V(min_y_ind,:) = [R*cos(min_alpha) R*sin(min_alpha) R]; % relocate the node to the origion.
% [~, min_y_ind]=min(abs(V(:,2)));
% V(min_y_ind,2) = 0;
% [~, min_z_ind]=min(V(:,3));
% V(min_z_ind,3) = 0;
%%
% Plotting meshed model
if plot_indenter
    cFigure; hold on;
    title('The meshed model','FontSize',fontSize);

    gpatch(F,V,'g');
    plotV(V1,'b-','LineWidth',2);
    plotV(V2,'b-','LineWidth',2);
    plotV(V3,'b-','LineWidth',2);

    axisGeom(gca,fontSize);
    camlight headlight;
    drawnow;
    axis equal;
end

% Make sure normals to tri3 surfaces of indentor point down towrads the
% specimen top surface
%     normal_2_surf = cross(V(F(end,2),:)-V(F(end,1),:),V(F(end,3),:)-V(F(end,1),:));
%     if(normal_2_surf(3)>1)
%         F = flip(F,2);
%     end
    
    
for i=1:size(F,1)
    V1 = V(F(i,1),:);
    V2 = V(F(i,2),:);
    V3 = V(F(i,3),:);
%     centroid = mean([V1;V2;V3]);
    normal = cross(V2-V1,V3-V1);
%     quiver3(centroid(1),centroid(2),centroid(3),normal(1),normal(2),normal(3));
    if normal(3)>0
        F(i,:) = flip(F(i,:),2);
    end
end
%%
%
% <<gibbVerySmall.gif>>
%
% _*GIBBON*_
% <www.gibboncode.org>
%
% _Kevin Mattheus Moerman_, <gibbon.toolbox@gmail.com>
% close all
% clear all
% clc
% 
% V = rand(3,3);
% P = [0.5 0.5 0.5]
% hold on
% scatter3(V(:,1),V(:,2),V(:,3))
% V = setCorner(V,P);
% scatter3(V(:,1),V(:,2),V(:,3), 'x')





function [V,I] = setCorner(V,P)
[~,I] = min((V(:,1)-P(1)).^2+(V(:,2)-P(2)).^2+(V(:,3)-P(3)).^2);
V(I,:) = P;
end

end

