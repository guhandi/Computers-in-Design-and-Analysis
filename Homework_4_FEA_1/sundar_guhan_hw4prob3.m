close all
clear all
clc

E = 210 * 10^9;
A = 25/10000;
P = 90000;
l1 = sqrt(2);
l2 = 1;
l3 = sqrt(2);
l4 = 1;
l5 = 1;

c1 = E*A/l1;
c2 = E*A/l2;
c3 = E*A/l3;
c4 = E*A/l4;
c5 = E*A/l5;

%rod1-2
b = 45;
c = cosd(b);
s = sind(b);

K1 = c1 * [c^2, c*s, -c^2, -c*s; c*s,s^2, -c*s, -s^2; -c^2, -c*s, c^2, c*s; -c*s, -s^2, c*s, s^2];

%rod2-4
b = 0;
c = cosd(b);
s = sind(b);

K2 = c2 * [c^2, c*s, -c^2, -c*s; c*s,s^2, -c*s, -s^2; -c^2, -c*s, c^2, c*s; -c*s, -s^2, c*s, s^2];

%rod4-3
b = 225;
c = cosd(b);
s = sind(b);

K3 = c3 * [c^2, c*s, -c^2, -c*s; c*s,s^2, -c*s, -s^2; -c^2, -c*s, c^2, c*s; -c*s, -s^2, c*s, s^2];

%rod3-1
b = 180;
c = cosd(b);
s = sind(b);

K4 = c4 * [c^2, c*s, -c^2, -c*s; c*s,s^2, -c*s, -s^2; -c^2, -c*s, c^2, c*s; -c*s, -s^2, c*s, s^2];

%rod3-2
b = 90;
c = cosd(b);
s = sind(b);

K5 = c5 * [c^2, c*s, -c^2, -c*s; c*s,s^2, -c*s, -s^2; -c^2, -c*s, c^2, c*s; -c*s, -s^2, c*s, s^2];

%Global Stiffness
KG = zeros(8,8);
for m=1:4
    K = zeros(8,8);
    switch m
        case 1
            Kn = K1;
            n1 = 1;
            n2 = 2;
        case 2
            Kn = K2;
            n1 = 2;
            n2 = 4;
        case 3
            Kn = K3;
            n1 = 4;
            n2 = 3;
        case 4
            Kn = K4;
            n1 = 3;
            n2 = 1;
        case 5
            Kn = K5;
            n1 = 3;
            n2 = 2;
        otherwise
            Kn = K1;
            end
    for i=1:4
        for j=1:4
       
            p = [2*n1-1, 2*n1, 2*n2-1, 2*n2];
            K(p(i),p(j)) =  Kn(i,j);
            KG = KG + K;
        end
    end
end

%forces
syms fx1 fx2 fx3 fx3 fx4 fy1 fy2 fy3 fy4
fx1 = 0;
fx2 = 0;
fx3 = 0;
fy2 = -90000;
fy3 = 0;

F = [fx1 fy1 fx2 fy2 fx3 fy3 fx4 fy4]';
KG_simp = KG([1,3,4,5,6], [1,3,4,5,6]);
F_simp = F([1,3,4,5,6]);
disp_simp = KG_simp\F_simp;

u1 = double(disp_simp(1));
u2 = double(disp_simp(2));
u3 = double(disp_simp(4));
u4=0;
v1 = 0;
v2 = double(disp_simp(3));
v3 = double(disp_simp(5));
v4 = 0;

displacement = [u1 v1 u2 v2 u3 v3 u4 v4]';
Force = KG * displacement;

%positions assuming node 1 is the origin
x(1) = 0;
y(1) = 0;
x(2) = 1;
y(2) = 1;
x(3) = 1;
y(3) = 0;
x(4) = 2;
y(4) = 1;

scale = 200;
for i=1:4
    x_deformed(i) = x(i) + scale * displacement(2*i-1);
    y_deformed(i) = y(i) + 200*displacement(2*i);
end

x_plot = [x(1) x(2) x(4) x(3) x(1) x(3) x(2)];
y_plot = [y(1) y(2) y(4) y(3) y(1) y(3) y(2)];
x_plotd = [x_deformed(1) x_deformed(2) x_deformed(4) x_deformed(3) x_deformed(1) x_deformed(3) x_deformed(2)];
y_plotd = [y_deformed(1) y_deformed(2) y_deformed(4) y_deformed(3) y_deformed(1) y_deformed(3) y_deformed(2)];

hold on
plot(x_plot, y_plot, 'LineWidth',2);
plot(x_plotd, y_plotd,'LineWidth',2);
xlabel('x (m)');
ylabel('y (m)');
title('Orignal vs "scaled" deformed truss structure');
legend('original', '"scaled" deformed');

disp('Global Stiffness');
disp(KG);

disp('Forces');
disp(Force);

disp('Displacements u1 v1 u2 v2...');
disp(displacement);


