function polar3D(fname,TSmin,TSmax)
% polar3D - 3D polar plot of bsl data with -60dB reference 
% reading 3 column ascii data with one header line
% MM 25.8.2025

%if(nargin<1) fname = '../out/mm-bem-bsl.txt'; end
%if(nargin<1) fname = '../out/sphere-1.905-510-bsl.txt'; end
%if(nargin<1) fname = '../out/soft-gypsilab-bsl.txt'; end
if(nargin<1) fname = '../out/YFT_swimbladder_origin-gypsilab-bsl.txt'; end
if(nargin<2) TSmin = -60; end
if(nargin<2) TSmax = -20; end

% reading file
f = fopen(fname);l=fgets(f);
d = fscanf(f,'%g');
n = length(d);
xyz = reshape(d, [3,n/3])';
fclose(f);

% transform to square matrices
n = sqrt(length(xyz));
theta = reshape(xyz(:,1)*pi/180, [n,n]);
phi = reshape(xyz(:,2)*pi/180, [n,n]);
TS = reshape(20*log10(xyz(:,3)), [n,n]);
TS(find(TS<TSmin)) = TSmin;
r = TS - TSmin;

% polar plot using xyz coordinates
X = cos(theta).*cos(phi).*r;
Y = cos(theta).*sin(phi).*r;
Z = sin(theta).*r;
surf(X,Y,Z,TS);
axis('equal');
caxis([TSmin,TSmax]);
ax = gca; ax.BoxStyle = 'full';
colorbar('Ticks',[TSmin:5:TSmax]);
xlabel('x');
ylabel('y');
zlabel('z');
hold on
p=40;
plot3([-p 0], [0 0], [0 0], 'r--', 'LineWidth', 1) % X-axis
plot3([0 p], [0 0], [0 0], 'r', 'LineWidth', 1) % X-axis
plot3([0 0], [-p 0], [0 0], 'g--', 'LineWidth', 1) % Y-axis
plot3([0 0], [0, p], [0 0], 'g', 'LineWidth', 1) % Y-axis
plot3([0 0], [0 0], [-p 0], 'b--', 'LineWidth', 1) % Z-axis
plot3([0 0], [0 0], [0 p], 'b', 'LineWidth', 1) % Z-axis
hold off
box off
%view(135,10);
shg

