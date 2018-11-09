%Problem 2


n = 17;
P1 = 10000;
P2 = 10000;
P3 = 10000;
E = 210 * 10^9;
A = 1.5 * 10^-3;


b = ones(1,n);

theta(1)= 60; 
theta(2)= 60; 
theta(3) = 120; 
theta(4) = 180; 
theta(5) = 120; 
theta(6) = 180;
theta(7) = 60; 
theta(8) = 180;
theta(9) = 120; 
theta(10)=180;
theta(11) = 60;
theta(12) = 180;
theta(13) = 120; 
theta(14) = 0; 
theta(15) = 60; 
theta(16) = 120; 
theta(17) = 60; 


for i=1:n
    c(i)=cosd(theta(i));
    s(i)=sind(theta(i));
end

syms S C
Kgeneral = [C^2 C*S -C^2 -C*S; C*S S^2 -C*S -S^2; -C^2 -C*S C^2 C*S;...
                -C*S -S^2 C*S S^2];
         
for i=1:n
    K{i} = double((A*E/b(i))*subs(Kgeneral, {C, S} , {c(i), s(i)})); 
end

 
%global K
N=20;
kg=zeros(N,N); 


for e=1:17
    for i=1:4
        for j=1:4
            switch e
                case 1
                    location = [1,2,3,4];
                case 2
                    location = [3,4,5,6];
                case 3
                    location = [2*3-1,2*3,2*4-1,2*4];
                case 4
                    location = [2*2-1,2*2,2*4-1,2*4];
                case 5
                    location = [2*2-1,2*2,2*5-1,2*5];
                case 6
                    location = [2*1-1,2*1,2*5-1,2*5];
                case 7
                    location = [2*4-1,2*4,2*5-1,2*5];
                case 8
                    location = [2*5-1,2*5,2*6-1,2*6];
                case 9
                    location = [2*4-1,2*4,2*6-1,2*6];
                case 10
                    location = [2*4-1,2*4,2*7-1,2*7];
                case 11
                    location = [2*6-1,2*6,2*7-1,2*7];
                case 12
                    location = [2*6-1,2*6,2*8-1,2*8];
                case 13
                    location = [2*7-1,2*7,2*8-1,2*8];
                case 14
                    location = [2*7-1,2*7,2*9-1,2*9];
                case 15
                    location = [2*7-1,2*7,2*10-1,2*10];
                case 16
                    location = [2*9-1,2*9,2*10-1,2*10];
                case 17
                    location = [2*8-1,2*8,2*9-1,2*9];
            end
            K_temp = zeros(N,N);
            K_temp(location(i),location(j)) = K{e}(i,j);
            kg = kg + K_temp;
                    
        end 
    end
end

 
Kg_r = kg([1, 2, 3, 4, 5, 6, 7, 8, 9, 11, 12, 13, 14, 17, 18, 19, 20],...
[1, 2, 3, 4, 5, 6, 7, 8, 9, 11, 12, 13, 14, 17, 18, 19, 20]); 

%Q1
disp('');
disp('Question 1')
disp('');
disp('Matrix 2-norm of global stiffness matrix:');
disp(norm(kg));
disp('');
disp('Matrix 2-norm of reduced stiffness matrix:');
disp(norm(Kg_r));
disp('');
disp('Condition # of global stiffness matrix:');
disp(cond(kg));
disp('');
disp('Condition # of reduced stiffness matrix:');
disp(cond(Kg_r));
disp('');



%forces
Fx = zeros(1,10);
Fy = zeros(1,10);
Fx(3) = -P2;
Fy(1) = -P1;
Fx(10) = P3; 

Forces = [Fx(1); Fy(1); Fx(2); Fy(2); Fx(3); Fy(3); Fx(4); Fy(4); Fx(5); Fy(5);...
    Fx(6); Fy(6); Fx(7); Fy(7); Fx(8); Fy(8); Fx(9); Fy(9); Fx(10); Fy(10)];

 
Force_r = Forces([1, 2, 3, 4, 5, 6, 7, 8, 9, 11, 12, 13, 14, 17, 18, 19, 20]); 
disp_r = Kg_r\Force_r; 

%displacements
u1 = double(disp_r(1));
v1 = double(disp_r(2)); 
u2 = double(disp_r(3));
v2 = double(disp_r(4));
u3 = double(disp_r(5));
v3 = double(disp_r(6)); 
u4= double(disp_r(7)); 
v4 = double(disp_r(8));
u5=double(disp_r(9));
v5 = 0; 
u6 = double(disp_r(10));
v6 = double(disp_r(11)); 
u7 = double(disp_r(12));
v7 = double(disp_r(13));
u8 = 0;
v8 = 0; 
u9 = double(disp_r(14));
v9 = double(disp_r(15));
u10 = double(disp_r(16));
v10 = double(disp_r(17)); 

displacement = [u1 v1 u2 v2 u3 v3 u4 v4 u5 v5 u6 v6 u7 v7 u8 v8 u9 v9 u10 v10]';
Force = kg * displacement; 

%Q2
disp(' ');
disp('Question 2');
max_disp = 0;
node_num_disp = 0;
max_stress = 0;
node_num_stress = 0;

for i=1:10
    
    x_disp = displacement(2*i-1);
    y_disp = displacement(2*i);
    total_disp = sqrt(x_disp^2 + y_disp^2);
    if (total_disp > max_disp)
        max_disp = total_disp;
        node_num_disp = i;
    end
    
    fx = Force(2*i-1);
    fy = Force(2*i);
    total_stress = sqrt(fx^2 + fy^2);
    if (total_stress > max_stress)
        max_stress = total_stress;
        node_num_stress = i;
    end
    
    
    disp(' ');
    label = ['Node ', num2str(i)];
    disp(label);
    Fx_label = ['External Fx: ', num2str(fx)];
    disp(Fx_label);
    Fy_label = ['External Fy: ', num2str(fy)];
    disp(Fy_label);

    x_label = ['Displacement in X: ', num2str(x_disp)];
    disp(x_label);
    y_label = ['Displacement in Y: ', num2str(y_disp)];
    disp(y_label);
end

%Q3
disp(' '); 
disp('Question 3');
disp(' ');
disp_lbl = ['Max displacement is: ', num2str(max_disp)];
disp(disp_lbl);
node_lbl = ['At node: ', num2str(node_num_disp)];
disp(node_lbl);

%Q4
disp(' '); 
disp('Question 4');
disp(' ');
stress_lbl = ['Max stress (tensile) is: ', num2str(max_stress)];
disp(stress_lbl);
node_lbl = ['At node: ', num2str(node_num_stress)];
disp(node_lbl);



%orginal coordinates
r3 = sqrt(3);

x(1) = -3;
y(1) = 0;
x(2) = -(2+1/2); 
y(2) = r3/2; 
x(3) = -2;
y(3) = r3; 
x(4) = -(1 + 1/2); 
y(4) = r3/2; 
x(5) = -2;
y(5) = 0;
x(6) = -1;
y(6) = 0; 
x(7) = -(1/2); 
y(7) = r3/2;
x(8) = 0;
y(8) = 0;
x(9) = 1/2;
y(9)= r3/2;
x(10)=0;
y(10)=r3; 


%deformed coordinates
scale = 200;
for i=1:10
    x_def(i) = x(i) + scale*displacement(2*i-1);   
    y_def(i) = y(i) + scale*displacement(2*i);  
end

original = [x(2)-x(1), x(3)-x(2), x(4)-x(3), x(4)-x(2), x(5)-x(2),x(5)-x(1), x(5)-x(4), x(6)-x(5), x(6)-x(4),x(7)-x(4),x(7)-x(6),x(8)-x(6),x(8)-x(7),x(9)-x(7),x(10)-x(7),x(10)-x(9),x(9)-x(8);...
          y(2)-y(1), y(3)-y(2), y(4)-y(3), y(4)-y(2), y(5)-y(2),y(5)-y(1), y(5)-y(4), y(6)-y(5), y(6)-y(4),y(7)-y(4),y(7)-y(6),y(8)-y(6),y(8)-y(7),y(9)-y(7),y(10)-y(7),y(10)-y(9),y(9)-y(8)];
     
deformed = [x_def(2)-x_def(1), x_def(3)-x_def(2), x_def(4)-x_def(3), x_def(4)-x_def(2), x_def(5)-x_def(2),x_def(5)-x_def(1), x_def(5)-x_def(4), x_def(6)-x_def(5), x_def(6)-x_def(4),x_def(7)-x_def(4),x_def(7)-x_def(6),x_def(8)-x_def(6),x_def(8)-x_def(7),x_def(9)-x_def(7),x_def(10)-x_def(7),x_def(10)-x_def(9),x_def(9)-x_def(8);...
    y_def(2)-y_def(1), y_def(3)-y_def(2), y_def(4)-y_def(3), y_def(4)-y_def(2), y_def(5)-y_def(2),y_def(5)-y_def(1), y_def(5)-y_def(4), y_def(6)-y_def(5), y_def(6)-y_def(4),y_def(7)-y_def(4),y_def(7)-y_def(6),y_def(8)-y_def(6),y_def(8)-y_def(7),y_def(9)-y_def(7),y_def(10)-y_def(7),y_def(10)-y_def(9),y_def(9)-y_def(8)];


for i=1:17
    strain(1,i)=(deformed(1,i)-original(1,i))/(original(1,i));
    strain(2,i)=(deformed(2,i)-original(2,i))/(original(2,i));
end

stress = E.*strain; 

x_plot = [x(1) x(2) x(1) x(5) x(2) x(3) x(2) x(4) x(2) x(5) x(4) x(3) x(4) x(5) x(4) x(7) x(4) x(6) x(5) x(6) x(6) x(7) x(6) x(8) x(7) x(8) x(7) x(9) x(7) x(10) x(8) x(9) x(10) x(9)];
y_plot = [y(1) y(2) y(1) y(5) y(2) y(3) y(2) y(4) y(2) y(5) y(4) y(3) y(4) y(5) y(4) y(7) y(4) y(6) y(5) y(6) y(6) y(7) y(6) y(8) y(7) y(8) y(7) y(9) y(7) y(10) y(8) y(9) y(10) y(9)];


% deformed
xdef_plot = [x_def(1) x_def(2) x_def(1) x_def(5) x_def(2) x_def(3) x_def(2) x_def(4) x_def(2) x_def(5) x_def(4) x_def(3) x_def(4) x_def(5) x_def(4) x_def(7) x_def(4) x_def(6) x_def(5) x_def(6) x_def(6) x_def(7) x_def(6) x_def(8) x_def(7) x_def(8) x_def(7) x_def(9) x_def(7) x_def(10) x_def(8) x_def(9) x_def(10) x_def(9)];
ydef_plot = [y_def(1) y_def(2) y_def(1) y_def(5) y_def(2) y_def(3) y_def(2) y_def(4) y_def(2) y_def(5) y_def(4) y_def(3) y_def(4) y_def(5) y_def(4) y_def(7) y_def(4) y_def(6) y_def(5) y_def(6) y_def(6) y_def(7) y_def(6) y_def(8) y_def(7) y_def(8) y_def(7) y_def(9) y_def(7) y_def(10) y_def(8) y_def(9) y_def(10) y_def(9)];

%plot
hold on
plot(x_plot, y_plot, 'LineWidth', 2);
plot(xdef_plot, ydef_plot, 'LineWidth', 2);

xlabel('x (meters)');
ylabel('y (meters)');
title('Original Truss vs. "Scaled" Deformed Truss');
legend('Original Truss', 'Deformed Truss (scaled 200x)');
hold off

