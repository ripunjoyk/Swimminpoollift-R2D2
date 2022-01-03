clc;
% clear variables;


Designation='50 X 25 X 2.9';%input('Section shape\n');

%% Mass

mg=0;
Mg=1200;
F_load=mg+Mg;
Fs=1500;
Sections=50;

%% Length Definition

Lac=.5905;
Lad=.493;
Lcd=.15207;
Lgh=.15207;
Lcg=.880;
Ldh=.880;
Lck=.665;
Lce=.60978;
Lgj=.30093;
delta=acosd(Lac/Lce);


%%

theta=linspace(70,140,100);
theta_2=(theta-90);
temp=theta-delta;
L=sqrt(Lce^2+Lck^2-2*Lce*Lck*cosd(temp));
phi=asind(Lce.*sind(temp)./L);
% EK=sqrt(Lck^2+Lce^2-2*Lck*Lce*cosd(temp));
% psi=asind(EK/(Lce*sind(temp)));

%% Forces

M_load=F_load*Lgj;
Fx_spring=Fs.*cosd(phi);
Fy_spring=Fs.*sind(phi);
Subtend=atand(110/105);
Fhx=(M_load./(Lgh.*sind(theta-Subtend)));
Fgx=+Fhx-F_load*cosd(180-theta);
Fgy=F_load*sind(180-theta);
Fdx=Fhx;
Fcx=-Fgx-Fx_spring;
Fcy=Fgy-Fy_spring;

%% Force and moment distribution

%For lower beam
F_axial_DH=-Fdx;
%For upper beam
%Axial
dL=Lcg/Sections;
l=dL;
for i=1:Sections
    if l<(Lcg-Lck)
        F_axial_CG(i,:)=Fgx;
    else
        F_axial_CG(i,:)=Fgx+Fx_spring;
    end
    l=l+dL;
end
%Shear
l=dL;
for i=1:Sections
    if l<(Lcg-Lck)
        F_shear_CG(i,:)=Fgy;
    else
        F_shear_CG(i,:)=Fgy-Fy_spring;
    end
    l=l+dL;
end
%Moment
l=dL;
for i=1:Sections
    if l<(Lcg-Lck)
        M_shear_CG(i,:)=Fgy*l;
    else
        M_shear_CG(i,:)=Fgy*l-Fy_spring*(Lck+l-Lcg);
    end
    l=l+dL;
end

F_shear_CG=flip(F_shear_CG);
M_shear_CG=flip(M_shear_CG);
F_axial_CG=flip(F_axial_CG);
l_array=linspace(0,Lcg,Sections);
%% Plot

figure(1)
hold on
plot(Fgx);
plot(Fgy);
% plot(Fhx);
% plot(Fdx);
plot(Fcx);
plot(Fcy);
plot(F_axial_DH);
legend('Fgx','Fgy','Fcx','Fcy','FaxialDH');
hold off

figure(2)
surf(theta,l_array,F_axial_CG);

figure(3)
surf(theta,l_array,F_shear_CG);

figure(4)
surf(theta,l_array,M_shear_CG);




%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
[Section,Ixx,SecArea,X,Y,Tx,Ty]=Data(Designation);
for i=1:100
    for j=1:Sections
        [Sigma_vm_max(j,i)]=Stress_P(F_shear_CG(j,i),F_axial_CG(j,i),M_shear_CG(j,i),0,Section,Ixx,SecArea,X,Y,Tx,Ty);
    end
end

A = max(Sigma_vm_max, [], 'all');

figure(5)
surf(theta,l_array,Sigma_vm_max);

%% 
function [Section,Ixx,SecArea,X,Y,Tx,Ty] = Data(Designation)
opts = spreadsheetImportOptions("NumVariables", 12);

% Specify sheet and range
opts.Sheet = "Sheet1";
opts.DataRange = "A2:L277";

% Specify column names and types
opts.VariableNames = ["Designation", "Section", "D", "B", "t", "T", "M", "A", "Ix", "Iy", "Rx", "Ry"];
opts.VariableTypes = ["string", "categorical", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double"];

% Specify variable properties
opts = setvaropts(opts, "Designation", "WhitespaceRule", "preserve");
opts = setvaropts(opts, ["Designation", "Section"], "EmptyFieldRule", "auto");

% Import the data
Crosssectionproperties = readtable("C:\Users\ripun\MTech_Project\data\Cross section properties.xlsx", opts, "UseExcel", false);
l=height(Crosssectionproperties);
clear opts
l=height(Crosssectionproperties);
for loop=1:l
    if Designation==Crosssectionproperties.Designation(loop)
        Section=Crosssectionproperties.Section(loop);
        Ixx=Crosssectionproperties.Ix(loop);
        SecArea=Crosssectionproperties.A(loop);
        X=Crosssectionproperties.B(loop);
        Y=Crosssectionproperties.D(loop);
        Tx=Crosssectionproperties.t(loop);
        Ty=Crosssectionproperties.T(loop);
    end
end

end
%% 
function [Sigma_vm_max] = Stress_P(F_shear_CG,F_a,M_shear_CG,M,Section,Ixx,SecArea,X,Y,Tx,Ty)
        Sigma_min = 0;
        Sigma_min_2 = 0;
%         [Section,Ixx,SecArea,X,Y,Tx,Ty]=Data(Designation);

switch Section
    case 'P'
        
        SecArea_1=SecArea*10^-4;
        Ixx_1=Ixx*10^-8;
        Sigma_T = F_a/SecArea_1;
        Sigma_M = M_shear_CG*(Y/2)/Ixx_1*10^(-3);
        Sigma_MM = M*(Y/2)/Ixx_1*10^(-3);
        Sigma_max = (Sigma_T+Sigma_M-Sigma_MM)/10^6;
        Sigma_max_2 = (Sigma_T-Sigma_M+Sigma_MM)/10^6;


        Sigma_max_In = Sigma_max;
        Sigma_min_In = 0;
        Sigma_max_In_2 = Sigma_max;
        Sigma_min_In_2 = 0;

    case 'B'
        
        SecArea_1=SecArea*10^-4;
        Ixx_1=Ixx*10^-8;
        Sigma_T = F_a/SecArea_1;
        Sigma_M = M_shear_CG*(Y/2)/Ixx_1*10^(-3);
        Sigma_MM = M*(Y/2)/Ixx_1*10^(-3);
        Sigma_max = (Sigma_T+Sigma_M-Sigma_MM)/10^6;
        Sigma_max_2 = (Sigma_T-Sigma_M+Sigma_MM)/10^6;
        
        Q = (Ty*X)*(Y/2-Ty/2)*10^(-9);
        tau = F_shear_CG*Q/Ixx_1/Tx/2*10^(3);
        Sigma_MM_2 = M*(Y/2-Ty)/Ixx_1*10^(-3);
        Sigma_M_temp = M_shear_CG*(Y/2-Ty)/Ixx_1*10^(-3)+Sigma_T-Sigma_MM_2;
        Sigma_M_temp_2 = -M_shear_CG*(Y/2-Ty)/Ixx_1*10^(-3)+Sigma_T+Sigma_MM_2;
        Sigma_max_In = (Sigma_M_temp/2+sqrt((Sigma_M_temp/2)^2+tau^2))/10^6;
        Sigma_min_In = (Sigma_M_temp/2-sqrt((Sigma_M_temp/2)^2+tau^2))/10^6;
        Sigma_max_In_2 = (Sigma_M_temp_2/2+sqrt((Sigma_M_temp_2/2)^2+tau^2))/10^6;
        Sigma_min_In_2 = (Sigma_M_temp_2/2-sqrt((Sigma_M_temp_2/2)^2+tau^2))/10^6;
    case 'L'
        
        SecArea_1=SecArea*10^-4;
        Ixx_1=Ixx*10^-8;
        Sigma_T = F_a/SecArea_1;
        Sigma_M = M_shear_CG*(Y/2)/Ixx_1*10^(-3);
        Y=Y-((X*Ty^2/Tx+Y^2-Ty^2)/(2*(Ty*X/Tx-Ty+Y)));
        Sigma_MM = M*(Y)/Ixx_1*10^(-3);
        Sigma_max = (Sigma_T+Sigma_M-Sigma_MM)/10^6;
        

        Sigma_max_In = Sigma_max;
        Sigma_min_In = 0;
%         Q = (Ty*Y)*(X/2-Ty/2)*10^(-9);
%         tau = F_p*Q/Ixx/Tx*10^(3);
%         Sigma_MM_2 = M*(X/2-Ty)/Ixx*10^(-3);
%         Sigma_M_temp = (F_p*L)*(X/2-Ty)/Ixx*10^(-3)+Sigma_T-Sigma_MM_2;
%         Sigma_max_2 = (Sigma_M_temp/2+sqrt((Sigma_M_temp/2)^2+tau^2))/10^6;
%         Sigma_min_2 = (Sigma_M_temp/2-sqrt((Sigma_M_temp/2)^2+tau^2))/10^6;


    case 'C'
       
        SecArea_1=SecArea*10^-4;
        Ixx_1=Ixx*10^-8;
        Sigma_T = F_a/SecArea_1;
        Sigma_M = M_shear_CG*(Y/2)/Ixx_1*10^(-3);
        Sigma_MM = M*(Y/2)/Ixx_1*10^(-3);
        Sigma_max = (Sigma_T+Sigma_M-Sigma_MM)/10^6;
        Sigma_max_2 = (Sigma_T-Sigma_M+Sigma_MM)/10^6;

        Q = (Ty*X)*(Y/2-Ty/2)*10^(-9);
        tau = F_shear_CG*Q/Ixx_1/Tx*10^(3);
        Sigma_MM_2 = M*(Y/2-Ty)/Ixx_1*10^(-3);
        Sigma_M_temp = M_shear_CG*(Y/2-Ty)/Ixx_1*10^(-3)+Sigma_T-Sigma_MM_2;
        Sigma_M_temp_2 = -M_shear_CG*(Y/2-Ty)/Ixx_1*10^(-3)+Sigma_T+Sigma_MM_2;
        Sigma_max_In = (Sigma_M_temp/2+sqrt((Sigma_M_temp/2)^2+tau^2))/10^6;
        Sigma_min_In = (Sigma_M_temp/2-sqrt((Sigma_M_temp/2)^2+tau^2))/10^6;
        Sigma_max_In_2 = (Sigma_M_temp_2/2+sqrt((Sigma_M_temp_2/2)^2+tau^2))/10^6;
        Sigma_min_In_2 = (Sigma_M_temp_2/2-sqrt((Sigma_M_temp_2/2)^2+tau^2))/10^6;

    case 'I'
       
        SecArea_1=SecArea*10^-4;
        Ixx_1 = Ixx*10^-8;
        Sigma_T = F_a/SecArea_1;
        Sigma_M = M_shear_CG*(Y/2)/Ixx_1*10^(-3);
        Sigma_MM = M*(Y/2)/Ixx_1*10^(-3);
        Sigma_max = (Sigma_T+Sigma_M-Sigma_MM)/10^6;
        Sigma_max_2 = (Sigma_T-Sigma_M+Sigma_MM)/10^6;

        Q = (Ty*X)*(Y/2-Ty/2)*10^(-9);
        tau = F_shear_CG*Q/Ixx_1/Tx*10^(3);
        Sigma_MM_2 = M*(Y/2-Ty)/Ixx_1*10^(-3);
        Sigma_M_temp = M_shear_CG*(Y/2-Ty)/Ixx_1*10^(-3)+Sigma_T-Sigma_MM_2;
        Sigma_M_temp_2 = -M_shear_CG*(Y/2-Ty)/Ixx_1*10^(-3)+Sigma_T+Sigma_MM_2;
        Sigma_max_In = (Sigma_M_temp/2+sqrt((Sigma_M_temp/2)^2+tau^2))/10^6;
        Sigma_min_In = (Sigma_M_temp/2-sqrt((Sigma_M_temp/2)^2+tau^2))/10^6;
        Sigma_max_In_2 = (Sigma_M_temp_2/2+sqrt((Sigma_M_temp_2/2)^2+tau^2))/10^6;
        Sigma_min_In_2 = (Sigma_M_temp_2/2-sqrt((Sigma_M_temp_2/2)^2+tau^2))/10^6;

end

Sigma_vm_1=sqrt(Sigma_max^2+Sigma_min^2-Sigma_max*Sigma_min);
Sigma_vm_2=sqrt(Sigma_max_2^2+Sigma_min_2^2-Sigma_max_2*Sigma_min_2);
Sigma_vm_3=sqrt(Sigma_max_In^2+Sigma_min_In^2-Sigma_max_In*Sigma_min_In);
Sigma_vm_4=sqrt(Sigma_max_In_2^2+Sigma_min_In_2^2-Sigma_max_In_2*Sigma_min_In_2);
Sigma_vm_max=max([Sigma_vm_1,Sigma_vm_2,Sigma_vm_3,Sigma_vm_4]);

end