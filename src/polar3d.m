function [theta,phi,r]=polar3D(fname,TSmin,TSmax)
% polar3D - 3D polar plot of bsl data with -60dB reference 
% reading 3 columns ascii data with one header line
% in first line there should be size pattern '# mxn'
% MM 25.8.2025

%fname = '../out/mm-bem-bsl.txt';
%fname = '../out/sphere-1.905-510-bsl.txt';
%fname = '../out/YFT_swimbladder_origin-gypsilab-bsl.txt';
%fname = '../out/YFT-mm-bem-bsl.txt';

%if(nargin<1) fname = '../out/soft-gypsilab-bsl.txt'; end
if(nargin<1) fname = '../out/YFT_swimbladder_origin-gypsilab-bsl.txt'; end
if(nargin<2) TSmin = -60; end
if(nargin<2) TSmax = -20; end

% reading file
f = fopen(fname);
line=fgets(f);
tokens = regexp(line, '#\s*(\d+)x(\d+)', 'tokens');
if ~isempty(tokens)
    m = str2double(tokens{1}{1});
    n = str2double(tokens{1}{2});
    %disp(sprintf('m = %d, n = %d\n', m, n));
else
    disp('Size pattern not found.'); return;
end
d = fscanf(f,'%g');
nn = length(d);
xyz = reshape(d, [3,nn/3])';
fclose(f);

% transform to square matrices
theta = reshape(xyz(:,1)*pi/180, [m,n]);
phi = reshape(xyz(:,2)*pi/180, [m,n]);
TS = reshape(20*log10(xyz(:,3)), [m,n]);
TS(find(TS<TSmin)) = TSmin;
r = TS - TSmin;

% polar plot using xyz coordinates
X = cos(theta).*cos(phi).*r;
Y = cos(theta).*sin(phi).*r;
Z = sin(theta).*r;
surf(Z,X,Y,TS);
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

