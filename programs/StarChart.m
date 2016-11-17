% filtering
stars = readtable('hyg3.csv');
vstars = stars(stars.mag <= 6.5, :);
cvstars = vstars(vstars.dist < 100000, :);
cvstars = sortrows(cvstars, 'mag');
cvstars = cvstars(2:end,:);

% inputs
JD = 2457754.5;
la = 38.6486;
lo = -90.3078;
r = 10;
ah = 2*pi/3;
av = 2*pi/3;

% compute azimuth A and altitude a
D = JD - 2451545.0;
GMST = 18.697374558 + 24.06570982441908*D;
L = 280.47 + 0.98565*D;
omega = 125.04 - 0.052954*D;
obliquity = 23.4393 - 0.0000004*D;
nutation = -0.000319*sind(omega) - 0.000024*sind(2*L);
GAST = GMST + nutation * cosd(obliquity);
LHA = (GAST - cvstars.ra) * 15 + lo;
a = asin(cosd(LHA).*cosd(cvstars.dec)*cosd(la) + sind(cvstars.dec)*sind(la));
A = atan2(-sind(LHA), (tand(cvstars.dec)*cosd(la)-sind(la)*cosd(LHA))) + pi;

% compute Cartesian coordinates
x = r*cos(a).*sin(A); x = x(a > 0);
y = r*cos(a).*cos(A); y = y(a > 0);
z = r*sin(a);         z = z(a > 0);
mag = cvstars.mag(a > 0);

% determine colors for 3D chart
c = linspace(1, .2, length(x));
c = [c' c' c'] + .2;

% 3D plotting
scatter3(x, y, z, 5, c, 'filled');
xlabel('X (East)');
ylabel('Y (North)');
zlabel('Z (Zenith)');
set(gca,'Color',[0 0 0]);
axis([-r,r,-r,r,0,r]);
title('3D Star Chart');

% compute projected 2D coordinates
% h*h + xl*xl + yl*yl = r*r;
h = r / sqrt(1 + tan(ah/2).^2 + tan(av/2).^2);
xl = h * tan(ah/2);
yl = h * tan(av/2);
x = x(z >= h); y = y(z >= h); mag = mag(z >= h); z = z(z >= h);
x = x*h./z; y = y*h./z;
y = y(x >= -xl & x <= xl); z = z(x >= -xl & x <= xl); mag = mag(x >= -xl & x <= xl); x = x(x >= -xl & x <= xl);
x = x(y >= -yl & y <= yl); z = z(y >= -yl & y <= yl); mag = mag(y >= -yl & y <= yl); y = y(y >= -yl & y <= yl);

% determine colors and dot sizes for 2D chart
mrange = max(mag) - min(mag);
mag = mag - min(mag);
sz = 5 + 10 * (1 - mag / mrange);
mag = 100.^.2.^(mag);
c = .8./mag+.2;
c = [c c c];

% 2D plotting
figure
scatter(x, y, sz, c, 'filled');
xlabel('E/(W)');
ylabel('N/(S)');
set(gca,'Color',[0 0 0]);
axis([-xl,xl,-yl,yl]);
title('2D Star Chart');