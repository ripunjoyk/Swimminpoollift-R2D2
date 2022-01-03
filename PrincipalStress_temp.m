clc;
clear variables;


L=0.1;%input('Length in metres\n');
F_a=-1000;%input('Axial force\n');
F_p=1000;%input('Perpendicular force\n');
M=-100;%input('Moment (+ for shagging)\n');
Section='25 X 25 X 2.6';%input('Section shape\n');
n=100;
dL=linspace(0.1,0,n);
for j=1:n


[Sigma_max,Sigma_min,Sigma_max_2,Sigma_min_2,Sigma_max_In,Sigma_min_In,Sigma_max_In_2,Sigma_min_In_2]=Stress_P(Section,F_a,F_p,M,dL(j));
Sigma_vm_max(j)=Stress_VM(Sigma_max,Sigma_min,Sigma_max_2,Sigma_min_2,Sigma_max_In,Sigma_min_In,Sigma_max_In_2,Sigma_min_In_2);

end

figure(1)
plot(dL,Sigma_vm_max);
xlabel('Length','FontSize',15);
ylabel('VM stress [N/mm^2]','FontSize',15);
legend('VM stress along the beam');

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
function [Sigma_max,Sigma_min,Sigma_max_2,Sigma_min_2,Sigma_max_In,Sigma_min_In,Sigma_max_In_2,Sigma_min_In_2] = Stress_P(Designation,F_a,F_p,M,L)
        Sigma_min = 0;
        Sigma_min_2 = 0;
        [Section,Ixx,SecArea,X,Y,Tx,Ty]=Data(Designation);

switch Section
    case 'P'
        
        SecArea=SecArea*10^-4;
        Ixx=Ixx*10^-8;
        Sigma_T = F_a/SecArea;
        Sigma_M = (F_p*L)*(Y/2)/Ixx*10^(-3);
        Sigma_MM = M*(Y/2)/Ixx*10^(-3);
        Sigma_max = (Sigma_T+Sigma_M-Sigma_MM)/10^6;
        Sigma_max_2 = (Sigma_T-Sigma_M+Sigma_MM)/10^6;


        Sigma_max_In = Sigma_max;
        Sigma_min_In = 0;
        Sigma_max_In_2 = Sigma_max;
        Sigma_min_In_2 = 0;

    case 'B'
        
        SecArea=SecArea*10^-4;
        Ixx=Ixx*10^-8;
        Sigma_T = F_a/SecArea;
        Sigma_M = (F_p*L)*(Y/2)/Ixx*10^(-3);
        Sigma_MM = M*(Y/2)/Ixx*10^(-3);
        Sigma_max = (Sigma_T+Sigma_M-Sigma_MM)/10^6;
        Sigma_max_2 = (Sigma_T-Sigma_M+Sigma_MM)/10^6;
        
        Q = (Ty*X)*(Y/2-Ty/2)*10^(-9);
        tau = F_p*Q/Ixx/Tx/2*10^(3);
        Sigma_MM_2 = M*(Y/2-Ty)/Ixx*10^(-3);
        Sigma_M_temp = (F_p*L)*(Y/2-Ty)/Ixx*10^(-3)+Sigma_T-Sigma_MM_2;
        Sigma_M_temp_2 = -(F_p*L)*(Y/2-Ty)/Ixx*10^(-3)+Sigma_T+Sigma_MM_2;
        Sigma_max_In = (Sigma_M_temp/2+sqrt((Sigma_M_temp/2)^2+tau^2))/10^6;
        Sigma_min_In = (Sigma_M_temp/2-sqrt((Sigma_M_temp/2)^2+tau^2))/10^6;
        Sigma_max_In_2 = (Sigma_M_temp_2/2+sqrt((Sigma_M_temp_2/2)^2+tau^2))/10^6;
        Sigma_min_In_2 = (Sigma_M_temp_2/2-sqrt((Sigma_M_temp_2/2)^2+tau^2))/10^6;
    case 'L'
        
        SecArea=SecArea*10^-4;
        Ixx=Ixx*10^-8;
        Sigma_T = F_a/SecArea;
        Sigma_M = (F_p*L)*(Y/2)/Ixx*10^(-3);
        Y=Y-((X*Ty^2/Tx+Y^2-Ty^2)/(2*(Ty*X/Tx-Ty+Y)));
        Sigma_MM = M*(Y)/Ixx*10^(-3);
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
       
        SecArea=SecArea*10^-4;
        Ixx=Ixx*10^-8;
        Sigma_T = F_a/SecArea;
        Sigma_M = (F_p*L)*(Y/2)/Ixx*10^(-3);
        Sigma_MM = M*(Y/2)/Ixx*10^(-3);
        Sigma_max = (Sigma_T+Sigma_M-Sigma_MM)/10^6;
        Sigma_max_2 = (Sigma_T-Sigma_M+Sigma_MM)/10^6;

        Q = (Ty*X)*(Y/2-Ty/2)*10^(-9);
        tau = F_p*Q/Ixx/Tx*10^(3);
        Sigma_MM_2 = M*(Y/2-Ty)/Ixx*10^(-3);
        Sigma_M_temp = (F_p*L)*(Y/2-Ty)/Ixx*10^(-3)+Sigma_T-Sigma_MM_2;
        Sigma_M_temp_2 = -(F_p*L)*(Y/2-Ty)/Ixx*10^(-3)+Sigma_T+Sigma_MM_2;
        Sigma_max_In = (Sigma_M_temp/2+sqrt((Sigma_M_temp/2)^2+tau^2))/10^6;
        Sigma_min_In = (Sigma_M_temp/2-sqrt((Sigma_M_temp/2)^2+tau^2))/10^6;
        Sigma_max_In_2 = (Sigma_M_temp_2/2+sqrt((Sigma_M_temp_2/2)^2+tau^2))/10^6;
        Sigma_min_In_2 = (Sigma_M_temp_2/2-sqrt((Sigma_M_temp_2/2)^2+tau^2))/10^6;

    case 'I'
       
        SecArea=SecArea*10^-4;
        Ixx = Ixx*10^-8;
        Sigma_T = F_a/SecArea;
        Sigma_M = (F_p*L)*(Y/2)/Ixx*10^(-3);
        Sigma_MM = M*(Y/2)/Ixx*10^(-3);
        Sigma_max = (Sigma_T+Sigma_M-Sigma_MM)/10^6;
        Sigma_max_2 = (Sigma_T-Sigma_M+Sigma_MM)/10^6;

        Q = (Ty*X)*(Y/2-Ty/2)*10^(-9);
        tau = F_p*Q/Ixx/Tx*10^(3);
        Sigma_MM_2 = M*(Y/2-Ty)/Ixx*10^(-3);
        Sigma_M_temp = (F_p*L)*(Y/2-Ty)/Ixx*10^(-3)+Sigma_T-Sigma_MM_2;
        Sigma_M_temp_2 = -(F_p*L)*(Y/2-Ty)/Ixx*10^(-3)+Sigma_T+Sigma_MM_2;
        Sigma_max_In = (Sigma_M_temp/2+sqrt((Sigma_M_temp/2)^2+tau^2))/10^6;
        Sigma_min_In = (Sigma_M_temp/2-sqrt((Sigma_M_temp/2)^2+tau^2))/10^6;
        Sigma_max_In_2 = (Sigma_M_temp_2/2+sqrt((Sigma_M_temp_2/2)^2+tau^2))/10^6;
        Sigma_min_In_2 = (Sigma_M_temp_2/2-sqrt((Sigma_M_temp_2/2)^2+tau^2))/10^6;

end
A=[Sigma_max,Sigma_min,Sigma_max_2,Sigma_min_2,Sigma_max_In,Sigma_min_In,Sigma_max_In_2,Sigma_min_In_2];
Temporary=max(A);
Temporary2=abs(min(A));
Temporary=max(Temporary,Temporary2);
disp(Temporary)
if Temporary==abs(Sigma_max)
        disp('Max stress occurs at the upper boundary and is given by :');
        disp(strcat(num2str(Sigma_max,'%.4f'),' MPa'));
elseif Temporary==abs(Sigma_max_2)
        disp('Max stress occurs at the lower boundary and is given by :');
        disp(strcat(num2str(Sigma_max_2,'%.4f'),' MPa'));
elseif Temporary==abs(Sigma_max_In)|| Temporary==abs(Sigma_min_In)
        disp('Max stress occurs at the point between the web and flange on the upper side and is given by');
        disp('Sigma_max');
        disp(strcat(num2str(Sigma_max_In,'%.4f'),' MPa'));
        disp('Sigma_min');
        disp(strcat(num2str(Sigma_min_In,'%.4f'),' MPa'));
elseif Temporary==abs(Sigma_max_In_2)|| Temporary==abs(Sigma_min_In_2)
        disp('Max stress occurs at the point between the web and flange on the lower side and is given by');
        disp('Sigma_max');
        disp(strcat(num2str(Sigma_max_In_2,'%.4f'),' MPa'));
        disp('Sigma_min');
        disp(strcat(num2str(Sigma_min_In_2,'%.4f'),' MPa'));
else
        disp('Something went wrong');
end
end
%%
function [Sigma_vm_max] = Stress_VM(Sigma_max,Sigma_min,Sigma_max_2,Sigma_min_2,Sigma_max_In,Sigma_min_In,Sigma_max_In_2,Sigma_min_In_2)

Sigma_vm_1=sqrt(Sigma_max^2+Sigma_min^2-Sigma_max*Sigma_min);
Sigma_vm_2=sqrt(Sigma_max_2^2+Sigma_min_2^2-Sigma_max_2*Sigma_min_2);
Sigma_vm_3=sqrt(Sigma_max_In^2+Sigma_min_In^2-Sigma_max_In*Sigma_min_In);
Sigma_vm_4=sqrt(Sigma_max_In_2^2+Sigma_min_In_2^2-Sigma_max_In_2*Sigma_min_In_2);
Sigma_vm_max=max([Sigma_vm_1,Sigma_vm_2,Sigma_vm_3,Sigma_vm_4]);

if Sigma_vm_max==abs(Sigma_vm_1)
        disp('Max VM stress occurs at the upper boundary and is given by :');

elseif Sigma_vm_max==abs(Sigma_vm_2)
        disp('Max VM stress occurs at the lower boundary and is given by :');

elseif Sigma_vm_max==abs(Sigma_vm_3)
        disp('Max VM stress occurs at the point between the web and flange on the upper side and is given by');

elseif Sigma_vm_max==abs(Sigma_vm_4)
        disp('Max VM stress occurs at the point between the web and flange on the lower side and is given by');

else
        disp('Something went wrong');
end
disp(strcat(num2str(Sigma_vm_max,'%.4f'),' MPa'));
end