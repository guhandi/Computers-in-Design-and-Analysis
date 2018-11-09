%Problem 1


L=4;
F_const=100;
Q=30;
M=2000;
E=200*10^9;

I(1) = (1/12)*(100*10^-3)*(100*10^-3)^3;
I(2) = (1/4)*pi*(30*10^-3)^4;
I(3) = (1/4)*pi*(20*10^-3)^4;

inertia = 5;

syms inertia len
K_general = [12*E*inertia/len^3  6*E*inertia/len^2 -12*E*inertia/len^3  6*E*inertia/len^2;...
    6*E*inertia/len^2  4*E*inertia/len  -6*E*inertia/len^2   2*E*inertia/len;...
    -12*E*inertia/len^3  -6*E*inertia/len^2  12*E*inertia/len^3  -6*E*inertia/len^2;...
    6*E*inertia/len^2  2*E*inertia/len  -6*E*inertia/len^2  4*E*inertia/len];


for j=1:16
    switch j
        case num2cell(1:4)
            I_act = I(1);
        case num2cell(5:12)
            I_act = I(2);
        case num2cell(13:16)
            I_act = I(3);
    end
               
    K{j} = subs(K_general,{inertia,len},{I_act,(L/4)});

end

%global k
Kg = zeros(34,34);
for i=1:16
    for j=1:4
        for k=1:4
            location=[2*i-1,2*i,2*i+1,2*i+2];
            Kg(location(j),location(k))=Kg(location(j),location(k))+ K{i}(j,k);
        end  
    end
end

%force and moments
F = zeros(1,34);
F(9)=Q*(L/4)/2;
F(10)=Q*(L/4)^2/12;
F(11)=Q*(L/4);
F(13)=Q*(L/4);
F(15)=Q*(L/4);
F(17)=Q*(L/4);
F(19)=Q*(L/4);
F(21)=Q*(L/4);
F(23)=Q*(L/4);
F(25)=Q*(L/4)/2-F_const;
F(26)=-Q*(L/4)^2/12;
F(34)=-M;

Kg_r=Kg([3:34],[3:34]);
F_r=F([3:34])';
disp_r = Kg_r\F_r;

u(1)=0;
u(2)=0;
for i=3:34
    u(i)=disp_r(i-2);
end  
Force_m = Kg*u';

%Q1
disp('Question 1')
disp(' ');
disp('Area Moment of Inertia');
disp('Element 1:');
disp(I(1));
disp('Element 2:');
disp(I(2));
disp('Element 3:');
disp(I(3));
disp(' ');


%Q2
disp('Question 2');
disp(' ');
disp('Matrix 2-norm of global stiffness matrix Kg:');
disp(norm(Kg));
disp(' ');
disp('Condition number of reduced stiffness matrix Kr:');
disp(cond(Kg_r));
disp(' ');

%Q3
disp('Question 3');
disp(' ');
for i=1:17
    disp(' ');
    node_lbl=['Node',num2str(i)];
    disp(node_lbl);
    disp_lbl = ['Displacement: ', num2str(u(2*i-1))];
    disp(disp_lbl);
    rot_lbl=['Rotation: ', num2str(u(2*i))];
    disp(rot_lbl);
end
disp(' ');

%Q4
disp('Question 4');
disp(' ');
rxnf_lbl = ['Reaction Force: ', num2str(Force_m(1))];
disp(rxnf_lbl);
rxnm_lbl=['Reaction Moment: ', num2str(Force_m(2))];
disp(rxnm_lbl);
disp(' ');

%Q5
disp('Question 5 & 6')
disp(' ');
figure(1);
for i=1:17
    v_deflect(i)=u(2*i-1);
end
length = [0:16];
plot(length,v_deflect);
xlabel('Length (meters)');
ylabel('Vertical Deflection (meters)');
title('Vertical Deflection vs Element Length');


figure(2);
for j=1:17
    rot_angle(j)=u(2*j)*(180/pi);
end
plot(length,rot_angle, 'k');
xlabel('Length (meters)');
ylabel('Rotation Angle (degrees)');
title('Rotation of Beam vs Element Length');