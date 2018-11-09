%1a
for i=1:5000
    P(i) = (100-2)+4*rand(1);
    L(i) = (0.5-0.01)+0.02*rand(1);
    E(i) = 200*10^9;
    I(i) = (0.1-0.005)+0.01*rand(1);
    S(i) = (P(i)*L(i)^3)/(E(i)*I(i));
end
figure
title('Uniform Distribution')
subplot(2,2,1)
histogram(P,60)
title('Load distribution')
xlabel('load (N)')
ylabel('# of samples')

subplot(2,2,2)
histogram(L,60)
title('Length of beam distribution')
xlabel('length (m)')
ylabel('# of samples')

subplot(2,2,3)
histogram(I,60)
title('Area moment of inertia distribution')
xlabel('area moment of inertia (m^4)')
ylabel('# of samples')

subplot(2,2,4)
histogram(S,60)
title('Beam deflection distribution')
xlabel('deflection (m)')
ylabel('# of samples')

%1b
for i=1:5000
    Pn(i) = 100+(2/3)*sqrt(2)*erfinv(2*rand(1)-1);
    Ln(i) = 0.5+(0.01/3)*sqrt(2)*erfinv(2*rand(1)-1);
    En(i) = 200*10^9;
    In(i) = 0.1+(0.005/3)*sqrt(2)*erfinv(2*rand(1)-1);;
    Sn(i) = (P(i)*L(i)^3)/(E(i)*I(i));
end

figure
title('Gaussian Distribution')
subplot(2,2,1)
histogram(Pn,60)
title('Load distribution')
xlabel('load (N)')
ylabel('# of samples')

subplot(2,2,2)
histogram(Ln,60)
title('Length of beam distribution')
xlabel('length (m)')
ylabel('# of samples')

subplot(2,2,3)
histogram(In,60)
title('Area moment of inertia distribution')
xlabel('area moment of inertia (m^4)')
ylabel('# of samples')

subplot(2,2,4)
histogram(Sn,60)
title('Beam deflection distribution')
xlabel('deflection (m)')
ylabel('# of samples')

n = [1:10:1000];
for i=1:length(n)
    sum_x2 = 0;
    sum_cosx = 0;
    for j=1:n(i)
        x(j) = rand(1);
        sum_x2 = sum_x2 + (x(j)^2);
        sum_cosx = sum_cosx + cos(pi*x(j));
    end
    est_x2(i) = sum_x2/n(i);
    est_cosx(i) = sum_cosx/n(i);
end
integral_x2 = 1/3;
integral_cosx = 0;


figure
title('Monte Carlo Integration')
subplot(2,1,1)
hold on
plot(n, est_x2);
plot(n,ones(length(n)) *integral_x2, 'LineWidth', 2);
title('f(x) = x^2')
xlabel('n')
ylabel('E(f(x)');
legend('Monte Carlo Estimator', 'Theoretical Integral')
hold off

subplot(2,1,2)
hold on
plot(n, est_cosx);
plot(n,ones(length(n)) *integral_cosx, 'LineWidth', 2);;
title('f(x) = cos(pi*x)')
xlabel('n')
ylabel('E(f(x)')
legend('Monte Carlo Estimator', 'Theoretical Integral')
hold off


%3a
theta  = [1:1:360];
w = 500*2*pi/60;
for i=1:60
    L=10;
    B=60;
    y(i) = 0.5*L*(1-cos(pi*i/B));
    v(i) = (pi*L*w/(2*B))*sin(pi*i/B);
    a(i) = (L/2) * ((pi*w/B)^2) * cos(pi*i/B);
end

for i=61:80
    y(i) = 10;
    v(i) = 0;
    a(i) = 0;
end

for i = 81:150
    L = 15;
    B = 70;
    y(i) = L*(((i-80)/B) - (1/(2*pi))*sin((2*pi*(i-80))/B)) + 10;
    v(i) = (L*w/B)*(1-cos((2*pi*(i-80))/B));
    a(i) = (2*pi*L) * ((w/B)^2) * sin((2*pi*(i-80))/B);
end

for i=151:200
    y(i) = 25;
    v(i) = 0;
    a(i) = 0;
end
 
for i = 201:300
    L = 25;
    B = 100;
    y(i) = 25 - L*(  (10*(i-200)^3)/(B^3) - (15*(i-200)^4)/(B^4) + (6*(i-200)^5)/(B^5));
    v(i) = L*(  (30*w*(i-200)^2)/(B^3) - (60*w*(i-200)^3)/(B^4) + (30*w*(i-200)^4)/(B^5));
    a(i) = L*(  (60*w^2*(i-200))/(B^3) - (180*w^2*(i-200)^2)/(B^4) + (120*w^2*(i-200)^3)/(B^5));
end

for i=301:360
    y(i) = 0;
    v(i) = 0;
    a(i) = 0;
end

hold on
plot(theta,y, 'LineWidth', 2);
plot(theta,v, 'LineWidth', 2);
plot(theta,a, 'LineWidth', 2);
legend('displacement (mm)', 'velocity (mm/s)','acceleration (mm/s^2)');
title('Cam Profiles');
xlabel('theta (degrees)');
ylabel('Profile magnitude');
hold off


%3b
a = 35;
for i=1:360
    p(i) = atan((v(i)/w)/(a+y(i)));
end
figure
plot(theta,p, 'LineWidth', 2);
xlabel('theta (degrees)');
ylabel('pressure angle (degrees)');
title('Pressure angle vs theta');

%3c
basex = 30*cosd(theta);
basey = 30*sind(theta);

primex = 35*cosd(theta);
primey = 35*sind(theta);

camx = (30+y).*cosd(theta);
camy = (30+y).*sind(theta);

pitchx = (35+y).*cosd(theta);
pitchy = (35+y).*sind(theta);

figure
hold on
plot(basex, basey, 'LineWidth', 2)
plot(primex, primey, 'LineWidth', 2)
plot(camx, camy, 'LineWidth', 2)
plot(pitchx, pitchy, 'LineWidth', 2)
xlabel('x (mm)')
ylabel('y (mm)')
title('Contour Curves')
legend('base circle', 'prime circle', 'cam contour', 'pitch curve')
hold off


%3d
for i=0:30:270
    rot = [cosd(i) sind(i) 0; -sind(i) cosd(i) 0; 0 0 1];
    for j=1:360
        cpoint = [camx(j); camy(j); 1];
        ppoint = [primex(j); primey(j); 1];
        
        rot_cam = rot * cpoint;
        rot_prime = rot * ppoint;
        
        rcamx(j) = rot_cam(1);
        rcamy(j) = rot_cam(2);
        
        rprimex(j) = rot_prime(1);
        rprimey(j) = rot_prime(2);
    end
    figure
    hold on
    plot(rcamx, rcamy);
    plot(rprimex, rprimey);
    xlabel('x (mm)')
    ylabel('y (mm)')
    str = int2str(i);
    txt = 'Contour orientation at theta = ';
    title(strcat(txt, str));
    hold off
end
    
        
        
        
   
    
