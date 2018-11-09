% 2D Newton Raphson Method
clear all
close all
clc
global L1 L2 L3 a b
global theta1
% Constants
L1= 7;
L2= 12;
L3= 15;
a= 18;
b= -5;
t1 = [0:360];

for i = 1:361
    theta1 = t1(i)*pi/180;
    %Initial guesses
    theta2 = 45*pi/180;
    theta3 = 80*pi/180;
    x0 = [theta2; theta3];
    error = 1;

    while error > 0.001
        f = constraints(x0);
        error = norm(f);
        D = jacobian(x0);
        delta_x = D\f;
        x0 = x0 - delta_x;
    end
    t2(i) = x0(1)*180/pi;
    t3(i) = x0(2)*180/pi;
 
end

%(1)
hold on
plot(t1,t3);
xlabel('input crank angle');
ylabel('19cm long link angle');
title('Crank vs Rocker link angle');
hold off


%2
Px = 5*cosd(t2+90)+7*cosd(t1);
Py = 5*sind(t2+90)+7*sind(t1);
figure
plot(Px, Py);
xlabel('x (cm)');
ylabel('y (cm)');
title('Path of Point P');

%3
Px2 = 5*cosd(t2+90);
Py2 = 5*sind(t2+90);
Qx = 4*cosd(t3)+12*cosd(t2);
Qy = 4*sind(t3)+12*sind(t2);

for i=1:361
    d(i) = sqrt((Qx(i)-Px2(i))^2+(Qy(i)-Py2(i))^2);
end

figure
plot(t1,d)
xlabel('crank angle');
ylabel('distance (cm)');
title('Distance betwenn P and Q');

%4
dt1 = 60*2*pi/60;
for i=1:361
    dt2(i) = dt1 * (L1/L2) * (sind(t3(i)-t1(i))/sind(t2(i)-t3(i)));
    dt3(i) = dt1 * (L1/L3) * (sind(t2(i)-t1(i))/sind(t2(i)-t3(i)));
end

figure
plot(t1,dt3)
xlabel('input crank angle');
ylabel('angular velocity (rad/s)');
title('Angular velocity of Rocker link');

%5
for i=1:361
   r = sqrt(Px(i)^2 + Py(i)^2); %Point P radius
   tv(i) = r*dt2(i); %translational velocity
end

figure
plot(t1,tv);
xlabel('input crank angle');
ylabel('velocity (cm/s)');
title('Translational velocity during rotation');



function f = constraints(x)
% Evaluate the equations
global L1 L2 L3 a b
global theta1
theta2 = x(1);
theta3 = x(2);
f1 = L1*cos(theta1) + L2*cos(theta2) - L3*cos(theta3) - a;
f2 = L1*sin(theta1) + L2*sin(theta2) - L3*sin(theta3) - b;
f = [f1 f2]';
end

function D = jacobian(x)
% Evaluate the Jacobian
global L1 L2 L3 a b
theta2 = x(1);
theta3 = x(2);
d11 = -L2*sin(theta2);
d12 = L3*sin(theta3);
d21 = L2*cos(theta2);
d22 = -L3*cos(theta3);
D = [d11 d12; d21 d22];
end
