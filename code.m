%ASSIGNMENT 4
%MUHAMMED RISWAN
%234103445

clc
clear all;
format long
%---------------------------------------------------------------------%
%---------------------------------------------------------------------%
%here change the value
fprintf("\n________________________________________________________________________________________________________________");
fprintf("\nEnter the mechanical propertes \n")
E1 =input('E1 in GPa (eg 38.6): ');;%gpa
E2 = input('E2 in GPa (eg 8.27): ');;%gpa
nu12 = input('NU12(eg 0.28): ');
G12 = input('G12 in GPa (eg 4.14): ');;%Gpa
fprintf("________________________________________________________________________________________________________________\n");
fprintf("\nENTER THE ULTIMATE STRENGTH\n")
sut1=input('SUT1 in MPa (eg 1062): ');;%mpa
suc1=input('SUC1 in MPa (eg 610): ');;%mpa
sut2=input('SUT2 in MPa (eg 31): ');;%mpa
suc2=input('SUC2 in MPa (eg 118): ');;%mpa
shearstress=input('SHEAR STRESS in MPa (eg 72): ');;%mpa
%N = [100; 0; 0; 0; 0; 0]; %nx only there like this enter [nx;ny;nxy;mx;my;mxy]
fprintf("________________________________________________________________________________________________________________\n");
NINPUT = input('\neg: 100,0,0,0,0,0 \nEnter like this enter  Nx,Ny,Nxy,Mx,My,Mxy: ','s');
NINPUT = strrep(NINPUT, ',', ' ');
N = sscanf(NINPUT, '%f');
fprintf("\n________________________________________________________________________________________________________________\n");
%material properties
% E1 = 38.6;%GPa
% E2 = 8.27;%GPa
% nu12 = 0.28;
% G12 = 4.14;%GPa
% %strength ultimate
% sut1=1062;%MPa
% suc1=610;%MPa
% sut2=31;%MPa
% suc2=118;%MPa
% shearstress=72;%MPa

% Laminate properties
% Angle of lamina
fprintf("\n________________________________________________________________________________________________________________");
ang = input('\neg: 0,45,-45,90,90,-45,45,0\nEnter the angles of lamina in degrees (theta_1, theta_2, theta_3... ):','s');
ang = strrep(ang, ',', ' ');
theta = sscanf(ang, '%f');
fprintf("\n________________________________________________________________________________________________________________\n");
%theta = [0,45,-45,90,90,-45,45,0]; % Fiber orientation angles in degrees for each layer
t =0.125; % Thickness of each layer mm
%---------------------------------------------------------------------------------------%
%---------------------------------------------------------------------------------------%
%upto here change accordingly


nu21 = (nu12*E2)/E1;
n = length(theta);
W=t*n;

% DIFFERENT HEIGHT Z
for i = 1:n+1
    Z(i) = -W/2 + (i-1) * t; 
end



 q11 = E1 / (1 - nu12 * nu21);
    q22 = E2 / (1 - nu12 * nu21);
    q12 = nu12 * E2 / (1 - nu12 * nu21);
    q66 = G12;
    q16=0;
    q26=0;
    q=[q11 q12 q16;q12 q22 q26;q16 q26 q66];

A=zeros(3,3);
B=zeros(3,3);
D=zeros(3,3);
for i=1:n %different laminar

    c = cosd(theta(i));
    s = sind(theta(i));
     Q11 = q11 * c ^ 4 + 2 * (q12 + 2 * q66) * s ^ 2 * c ^ 2 + q22 * s ^ 4;
    Q12 = (q11 + q22 - 4 * q66) * c ^ 2 * s ^ 2 + q12 * (s ^ 4 + c ^ 4);
    Q22 = q11 * s ^ 4 + 2 * (q12 + 2 * q66) * c ^ 2 * s ^ 2 + q22 * c ^ 4;
    Q16 = (q11 - q12 - 2 * q66) * s * c ^ 3 - (q22 - q12 - 2 * q66) * s ^ 3 * c;
    Q26 = (q11 - q12 - 2 * q66) * s ^ 3 * c - (q22 - q12 - 2 * q66) * s * c ^ 3;
    Q66 = (q11 + q22 - 2 * q12 - 2 * q66) * c ^ 2 * s ^ 2 + q66 * (s ^ 4 + c ^ 4);
     Q_bars = [Q11 Q12 Q16; Q12 Q22 Q26;Q16 Q26 Q66];
      A=A+Q_bars*(Z(i+1)-Z(i));
      
B=B+Q_bars*(Z(i+1)^2-Z(i)^2)*0.5;
D=D+Q_bars*(Z(i+1)^3-Z(i)^3)/3;
end
invER=inv(A);
Ex=1/(W*invER(1,1));
Ey=1/(W*invER(2,2));
txy=-invER(1,2)/invER(1,1);
fprintf("\nENGINEERING CONSTANT EX=%f\n",Ex);
fprintf("ENGINEERING CONSTANT EY=%f\n",Ey);
fprintf("ENGINEERING CONSTANT TXY=%f\n\n",txy);
ABBD=[A B;B D];
 fprintf("ABBD matrix");
 ABBD
% Calculate the MID strain matrix 
Eo = inv(ABBD) * N;
 
E = Eo(1:3); % Reference strain mid
K = Eo(4:6); % Curvature mid
sigma1=zeros(3,n);
for i=1:n %strain
Ex= (E + Z(i)*K);%strain local
c = cosd(theta(i));
    s = sind(theta(i));
 T=[c^2 s^2 2*s*c;s^2 c^2 -2*s*c;-s*c s*c c^2-s^2];
E1=T*Ex;
E1=[E1(1);E1(2);2*E1(3)];
sigma1(:,i)=q*E1;
end


ratio=zeros(3,n);
for i=1:n
    if sigma1(1,i) > 0
        ratio(1,i)=sigma1(1,i)/sut1;
    else 
         ratio(1,i)=sigma1(1,i)/suc1;
    end
    if sigma1(2,i) > 0
        ratio(2,i)=sigma1(2,i)/sut2;
    else 
         ratio(2,i)=sigma1(2,i)/suc2;
    end
    if sigma1(3,i)>= 0
        ratio(3,i)=sigma1(3,i)/shearstress;
    end
end


fprintf("_____________________________________________________________________________________\n")
fprintf("lamina\t\tsigma1\t\tStress ratio\tsigma2\tStress ratio\tshear\tStress ratio\n");
fprintf("_____________________________________________________________________________________")
for i=1:n
    fprintf("\n\t%d\t\t%f\t%f\t\t%f\t%f\t%f\t%f\n",i,sigma1(1,i),ratio(1,i),sigma1(2,i),ratio(2,i),sigma1(3,i),ratio(3,i));
end
fprintf("_____________________________________________________________________________________\n")





% Finding the First ply failure mode
[~, idx] = max(abs(ratio(:)));
[sigma, lamina] = ind2sub(size(ratio), idx);

sigma1s=sigma1(sigma, lamina);
if sigma1s > 0 && sigma==1
        fprintf("\nfailure due to longitudinal tensile in %d degree lamina",theta(lamina));
    elseif sigma1s < 0 && sigma==1
        fprintf("\nfailure due to longitudinal compression in %d degree lamina",theta(lamina));
    elseif sigma1s > 0 && sigma==2
        fprintf("\nfailure due to transverse tensile in %d degree lamina",theta(lamina));
    elseif sigma1s < 0 && sigma==2
        fprintf("\nfailure due to transverse compression in %d degree lamina",theta(lamina));
    elseif sigma1s > 0  && sigma1s < 0 && sigma==3
       fprintf("\nfailure due to shear strength in %d degree lamina",theta(lamina));
    end


N_x = N(1) / ratio(sigma, lamina);
fprintf('\n\nif there is only load is there then:- \n%d degree layer will fail with FPF load = %f N/mm\n',theta(lamina),N_x);


fprintf("\n_____________________________________________________________\n");
% Now if laminate is subjected to temp change
fprintf("\ntemperature is there then enter temperature\n");
Tdiff=input('Temp difference in celcius (eg 50): ');;%mpa
%temperature is there then change this value
alpha1=input('Enter alpha_1 values in *10^-6  m/m/C  (eg 8.6): ');
alpha2=input('Enter alpha_2 values in *10^-6  m/m/C  (eg 22.1): ');
alpha12=input('Enter alpha_12 values in *10^-6  m/m/C  (eg 0): ');
alpha_1=alpha1*1.0e-6;%m/m/c
alpha_2=alpha2*1.0e-6;%m/m/c
alpha_12=alpha12*1.0e-6;%m/m/c
fprintf("\n_____________________________________________________________\n");

alpha=[alpha_1;alpha_2;alpha_12];
NXT=zeros(3,1);
MXT=zeros(3,1);
for i=1:n
    c = cosd(theta(i));
    s = sind(theta(i));
 T=[c^2 s^2 2*s*c;s^2 c^2 -2*s*c;-s*c s*c c^2-s^2];
alpha_x= inv(T)*alpha;
alpha_x=[alpha_x(1);alpha_x(2);2*alpha_x(3)];
 Q11 = q11 * c ^ 4 + 2 * (q12 + 2 * q66) * s ^ 2 * c ^ 2 + q22 * s ^ 4;
    Q12 = (q11 + q22 - 4 * q66) * c ^ 2 * s ^ 2 + q12 * (s ^ 4 + c ^ 4);
    Q22 = q11 * s ^ 4 + 2 * (q12 + 2 * q66) * c ^ 2 * s ^ 2 + q22 * c ^ 4;
    Q16 = (q11 - q12 - 2 * q66) * s * c ^ 3 - (q22 - q12 - 2 * q66) * s ^ 3 * c;
    Q26 = (q11 - q12 - 2 * q66) * s ^ 3 * c - (q22 - q12 - 2 * q66) * s * c ^ 3;
    Q66 = (q11 + q22 - 2 * q12 - 2 * q66) * c ^ 2 * s ^ 2 + q66 * (s ^ 4 + c ^ 4);
     Q_bars = [Q11 Q12 Q16; Q12 Q22 Q26;Q16 Q26 Q66];
      NXT=NXT+Tdiff*Q_bars*(Z(i+1)-Z(i))*alpha_x;
MXT=MXT+Tdiff*0.5*Q_bars*(Z(i+1)^2-Z(i)^2)*alpha_x;
end

NXTT=[NXT;MXT];
% Calculate the MID strain matrix  OF TEMPERATURE
EoT = inv(ABBD) * NXTT;
% Extract E and K from strain matrix e
ET = EoT(1:3); % Reference strain mid
KT = EoT(4:6); % Curvature mid
sigma1r=zeros(3,n);
strainr=zeros(3,n);
for i=1:n %strain
Ex= (ET + Z(i)*KT);%strain local
  c = cosd(theta(i));
    s = sind(theta(i));
 T=[c^2 s^2 2*s*c;s^2 c^2 -2*s*c;-s*c s*c c^2-s^2];
alpha_x= inv(T)*alpha;
alpha_x=[alpha_x(1);alpha_x(2);2*alpha_x(3)];
Ext=Tdiff*alpha_x;
Exr=Ex-Ext;
strainr(:,i)=Exr;
 Q11 = q11 * c ^ 4 + 2 * (q12 + 2 * q66) * s ^ 2 * c ^ 2 + q22 * s ^ 4;
    Q12 = (q11 + q22 - 4 * q66) * c ^ 2 * s ^ 2 + q12 * (s ^ 4 + c ^ 4);
    Q22 = q11 * s ^ 4 + 2 * (q12 + 2 * q66) * c ^ 2 * s ^ 2 + q22 * c ^ 4;
    Q16 = (q11 - q12 - 2 * q66) * s * c ^ 3 - (q22 - q12 - 2 * q66) * s ^ 3 * c;
    Q26 = (q11 - q12 - 2 * q66) * s ^ 3 * c - (q22 - q12 - 2 * q66) * s * c ^ 3;
    Q66 = (q11 + q22 - 2 * q12 - 2 * q66) * c ^ 2 * s ^ 2 + q66 * (s ^ 4 + c ^ 4);
     Q_bars = [Q11 Q12 Q16; Q12 Q22 Q26;Q16 Q26 Q66];
sigmar=Q_bars*Exr;
sigma1r(:,i)=T*sigmar*10^3;%here give
end
fprintf("\nRESIDUAL STRAINS IN PLIES IN XY CORDINATES\n");
fprintf("______________________________________________________________________\n");
fprintf("lamina\t\tstrain1 residual\t\t strain2 residual\t\tshear residual\t\n");
fprintf("______________________________________________________________________");
for i=1:n
    fprintf("\n\t%d\t\t\t%f\t\t\t%f\t\t\t%f\n",i,strainr(1,i),strainr(2,i),strainr(3,i));
end
fprintf("______________________________________________________________________\n")
fprintf("\nRESIDUAL STRESS IN PLIES IN 1-2 CORDINATES\n");
fprintf("______________________________________________________________________\n");
fprintf("lamina\t\tsigma1 residual\t\t sigma2 residual\t\tshear residual\t\n");
fprintf("______________________________________________________________________");
for i=1:n
    fprintf("\n\t%d\t\t\t%f\t\t\t%f\t\t\t%f\n",i,sigma1r(1,i),sigma1r(2,i),sigma1r(3,i));
end
fprintf("______________________________________________________________________\n")
sigma1re=sigma1r(sigma, lamina);
sigma1s=sigma1(sigma, lamina) ;
    fprintf('\nFor Temp. change of %d degree celcius, Residual thermal stress in %d degree layer: %f MPa\n\n',Tdiff,theta(lamina),sigma1re );
  % Finding the total stress in concerned ply
    Sigma = sigma1s + sigma1re;
   % Calculation for FPF load for thermal load
    if sigma1s > 0 && sigma==1
        su = sut1;
    elseif sigma1s < 0 && sigma==1
        su = suc1;
    elseif sigma1s > 0 && sigma==2
        su = sut2;
    elseif sigma1s < 0 && sigma==2
        su = suc2;
    elseif sigma1s > 0  && sigma1s < 0 && sigma==3
        su = shearstress;
    end
      sigma2 = su - (sigma1re);
    % FPF load at which sigma2 will attain its value
    Nxt = N(1)/(sigma1s/sigma2);

    fprintf('\n%d degree layer(s) will fail with temperature change of %d at FPF load = %f N/mm\n',theta(lamina),Tdiff,Nxt);

