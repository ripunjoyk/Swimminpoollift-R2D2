clc;
% clear variables;
clear all;

Designation='50 X 25 X 2.9';%input('Section shape\n');

%% Mass

mg=0;
Mg=1200;
F_load=mg+Mg;
Fs=1500;
Sections=150;
Angle_N=100;
%% Length Definition

Lac=.5905;
Lad=.493;
Lcdx=.15207;
Lgh=.15207;
Lcg=.880;
Ldh=.880;
Lck=.665;
Lce=.60978;
Lgj=.30093;
delta=acosd(Lac/Lce);
Lcdy=Lac-Lad;


%% MISCELLANEOUS

Strength=540; %% For mild steel
FoS=5; %% Assumed
Expected_Strength=Strength/FoS;
theta=linspace(70,140,Angle_N);
theta_2=(theta-90);
temp=theta-delta;
L=sqrt(Lce^2+Lck^2-2*Lce*Lck*cosd(temp));
phi=asind(Lce.*sind(temp)./L);

%% Forces

M_load=F_load*Lgj;
Fx_spring=Fs.*cosd(phi);
Fy_spring=Fs.*sind(phi);
Subtend=atand(110/105);
Fhx=(M_load./(Lgh.*sind(theta-Subtend)));
Fgx=+Fhx-F_load*cosd(180-theta);
Fgy=F_load*sind(180-theta);
Fdx=Fhx;
Fcx=Fgx+Fx_spring;
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
[Designation_temp,Section,Ixx,SecArea,X,Y,Tx,Ty]=Data(Designation,0,0);

dl=Lcg/Sections;
for i=1:Angle_N
    L=0;
    for j=1:Sections
        if L==0
            switch_hole=1;
            D=50/2;
        elseif abs(L-Lck)<=Lcg/(Sections-20)
            switch_hole=1;
            D=50/2;
        else
            D=0;
            switch_hole=0;
        end
        [Sigma_vm_max(j,i)]=Stress_P(F_shear_CG(j,i),F_axial_CG(j,i),M_shear_CG(j,i),0,Section,Ixx,SecArea,X,Y,Tx,Ty,switch_hole,D);
        L=L+dl;
    end
end

A = max(Sigma_vm_max, [], 'all');

figure(5)
surf(theta,l_array,Sigma_vm_max);

%%
%Axial
for i=1:Sections
        F_axial_AC(i,:)=Fcx.*sind(theta_2)-Fcy.*cosd(theta_2);
end

%Shear
for i=1:Sections
        F_shear_AC(i,:)=Fcy.*sind(theta_2)+Fcx.*cosd(theta_2);
end

%Moment
dL=Lad/Sections;
l=dL;
for i=1:Sections
    M_shear_AC(i,:)=Fcy.*cosd(theta_2)*l+Fcx.*sind(theta_2)*l;
    l=l+dL;
end
%Pure Moment
L_diag_base=sqrt(Lcdx^2+Lcdy^2);
theta_temp=acosd(Lcdy/L_diag_base);
temp_angle=90-theta_temp+theta_2;
for i=1:Sections
    M_AC(i,:)=-Fcy*L_diag_base.*cosd(temp_angle)+Fcx*L_diag_base.*sind(temp_angle);
end

F_shear_AC=flip(F_shear_AC);
M_shear_AC=flip(M_shear_AC);
F_axial_AC=flip(F_axial_AC);
M_AC=flip(M_AC);

D=0;%% Don't change
switch_hole=0;%% Don't change

for k=1:83
    [Designation_temp,Section,Ixx,SecArea,X,Y,Tx,Ty]=Data(Designation,1,k);
    for i=1:Angle_N
        for j=1:Sections
            [Sigma_vm_max(j,i)]=Stress_P(F_shear_AC(j,i),F_axial_AC(j,i),M_shear_AC(j,i),M_AC(j,i),Section,Ixx,SecArea,X,Y,Tx,Ty,switch_hole,D);
        end
    end
    A = max(Sigma_vm_max, [], 'all');
    if A<=Expected_Strength
        X=[Designation_temp,A];
        disp(X);
    end
end



  

%% 
function [Designation,Section,Ixx,SecArea,X,Y,Tx,Ty] = Data(Designation,Temp_switch,loop)
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
if Temp_switch==0
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
else
         Designation=Crosssectionproperties.Designation(loop+140);
         Section=Crosssectionproperties.Section(loop+140);
         Ixx=Crosssectionproperties.Ix(loop+140);
         SecArea=Crosssectionproperties.A(loop+140);
         X=Crosssectionproperties.B(loop+140);
         Y=Crosssectionproperties.D(loop+140);
         Tx=Crosssectionproperties.t(loop+140);
         Ty=Crosssectionproperties.T(loop+140);
end

end
%% 
function [Sigma_vm_max] = Stress_P(F_shear_CG,F_a,M_shear_CG,M,Section,Ixx,SecArea,X,Y,Tx,Ty,switch_hole,D)
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
        if switch_hole==1
            SecArea_1=SecArea-2*Tx*2*D/100;
            Ixx_1=Ixx-2*Tx*D^3/12;
            F_a=2.8*F_a;
            F_shear_CG=2.8*F_shear_CG;
            M_shear_CG=M_shear_CG*2*Y/D;
        else
            SecArea_1=SecArea;
            Ixx_1=Ixx;
        end
        SecArea_1=SecArea_1*10^-4;
        Ixx_1=Ixx_1*10^-8;
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