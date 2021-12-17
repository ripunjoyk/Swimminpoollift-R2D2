clc
clear all
close all
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

clear opts




[A,B,C] = crossectional_analysis('P',127,127,5.4,5.4,Crosssectionproperties);
disp('Pipe of dia P15')
disp(strcat('A = ',num2str(C,'%.2f')))
disp(strcat('Ixx = ',num2str(A,'%.2f')))
[A,B,C] = crossectional_analysis('B',60,40,2.9,2.9,Crosssectionproperties);
disp('Box of 60 x 40')
disp(strcat('A = ',num2str(C,'%.2f')))
disp(strcat('Ixx = ',num2str(A,'%.2f')))
[A,B,C] = crossectional_analysis('C',75,40,4.8,7.5,Crosssectionproperties);
disp('C channel of 75 x 40')
disp(strcat('A = ',num2str(C,'%.2f')))
disp(strcat('Ixx = ',num2str(A,'%.2f')))
[A,B,C] = crossectional_analysis('I',100,100,6,10,Crosssectionproperties);
disp('I channel of 75 x 45')
disp(strcat('A = ',num2str(C,'%.2f')))
disp(strcat('Ixx = ',num2str(A,'%.2f')))

function [Ixx,R,SecArea] = crossectional_analysis(Section,X,Y,Tx,Ty,Crosssectionproperties)
switch Section
    case 'P'
        SecArea = (pi*X*Y-pi*(X-2*Tx)*(Y-2*Ty))/400;
        Ixx = pi*(X*Y*X*Y-(X-2*Tx)^4)/640000;
        R = 20*Ixx/X;

    case 'B'
        SecArea = 2*Tx*((X-4*Tx)+(Y-4*Tx)+(3*pi*Tx/2))/100;
        Ixx =  ((X*Y*Y*Y)-(X-2*Tx)*(Y-2*Ty)^3)/120000;
        %         Ixx =  (T*(Y-4*T)^3/6)+ (( ((B-4T)*T^3/3)+(T*(B-4T)*(D-T)^2))/2) + (pi*t^4(405-(3136/pi^2))/108)+(3*pi*T^2)*((9*pi*(Y-4*T)+56*T)/18*pi)^2;
        R = 20*Ixx/Y;

    case 'L'
        SecArea = ((X-Tx)*Ty+(Y-Ty)*Tx+Tx*Ty)/100;
        Ixx = ((Tx*(Y-Ty)^3/12)+(X*Ty^3/12)+(X*Ty*(X-Ty)*(X-Ty)/4))/10000;
        R = 20*Ixx/Y;

    case 'C'
        SecArea = (2*X*Ty+(Y-2*Ty)*Tx)/100;
        Ixx = ((Tx*(Y-2*Ty)^3/12)+2*((Ty*Ty*Ty*X/12)+(Ty*X*(Y-Ty)*(Y-Ty)/4)))/10000;
        R = 20*Ixx/Y;

    case 'I'
        SecArea = (2*X*Ty+(Y-2*Ty)*Tx)/100;
        Ixx = (Tx*(Y-2*Ty)^3/12)+2*((Ty*Ty*Ty*X/12)+(Ty*X*(Y-Ty)*(Y-Ty)/4));
        Ixx = Ixx/10000;
        R = 20*Ixx/Y;

end

A = find(Crosssectionproperties.Section==Section & Crosssectionproperties.A<=SecArea*1.15 & Crosssectionproperties.A>=SecArea/1.15);
disp(Crosssectionproperties(A,:));
end
