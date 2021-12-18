clc;
%FS=input('Factor of safety\n');
L=input('Length in metres\n');
F_a=input('Axial force\n');
F_p=input('Perpendicular force\n');
Section=input('Section shape\n');
X=input('D\n');
Y=input('B\n');
Tx=input('t\n');
Ty=input('T\n');
[Sigma_max,Sigma_min,Sigma_max_2,Sigma_min_2]=Stress_P(Section,X,Y,Tx,Ty,F_a,F_p,L);
if Sigma_max>=Sigma_max_2
    disp('Max stress occurs at the boundary and is given by :');
    disp(strcat(num2str(Sigma_max,'%.4f'),' MPa'));
else
    disp('Max stress occurs at the point between the web and flange and is given by');
    disp('Sigma_max');
    disp(strcat(num2str(Sigma_max_2,'%.4f'),' MPa'));
    disp('Sigma_min');
    disp(strcat(num2str(Sigma_min_2,'%.4f'),' MPa'));
end

%%

function [Sigma_max,Sigma_min,Sigma_max_2,Sigma_min_2] = Stress_P(Section,X,Y,Tx,Ty,F_a,F_p,L)
switch Section
    case 'P'
        SecArea = (pi*X*Y-pi*(X-2*Tx)*(Y-2*Ty))/400*10^(-4);
        Ixx = pi*(X*Y*X*Y-(X-2*Tx)^4)/640000*10^(-8);
%         ElstMod = 20*Ixx/X;
        Sigma_T = F_a/SecArea;
        Sigma_M = (F_p*L)*(X/2)/Ixx*10^(-3);
        Sigma_max = (Sigma_T+Sigma_M)/10^6;
        Sigma_min = 0;

        Sigma_max_2 = Sigma_max;
        Sigma_min_2 = 0;        

    case 'B'
        SecArea = 2*Tx*((X-4*Tx)+(Y-4*Tx)+(3*pi*Tx/2))/100*10^(-4);
        Ixx =  ((X*Y*Y*Y)-(X-2*Tx)*(Y-2*Ty)^3)/120000*10^(-8);
%         ElstMod = 20*Ixx/Y;
        Sigma_T = F_a/SecArea;
        Sigma_M = (F_p*L)*(X/2)/Ixx*10^(-3);
        Sigma_max = (Sigma_T+Sigma_M)/10^6;
        Sigma_min = 0;
        
        Q = (Ty*Y)*(X/2-Ty/2)*10^(-9);
        tau = F_p*Q/Ixx/Tx/2*10^(3);
        Sigma_M_temp = (F_p*L)*(X/2-Ty)/(210*10^(9)*Ixx)*10^(-3)+Sigma_T;
        Sigma_max_2 = (Sigma_M_temp/2+sqrt((Sigma_M_temp/2)^2+tau^2))/10^6;
        Sigma_min_2 = (Sigma_M_temp/2-sqrt((Sigma_M_temp/2)^2+tau^2))/10^6;
    case 'L'
        SecArea = ((X-Tx)*Ty+(Y-Ty)*Tx+Tx*Ty)/100*10^(-4);
        Ixx = ((Tx*(Y-Ty)^3/12)+(X*Ty^3/12)+(X*Ty*(X-Ty)*(X-Ty)/4))/10000*10^(-8);
%         ElstMod = 20*Ixx/Y;
        Sigma_T = F_a/SecArea;
        Sigma_M = (F_p*L)*(X/2)/Ixx*10^(-3);
        Sigma_max = (Sigma_T+Sigma_M)/10^6;
        Sigma_min = 0;

        Q = (Ty*Y)*(X/2-Ty/2)*10^(-9);
        tau = F_p*Q/Ixx/Tx*10^(3);
        Sigma_M_temp = (F_p*L)*(X/2-Ty)/(210*10^(9)*Ixx)*10^(-3)+Sigma_T;
        Sigma_max_2 = (Sigma_M_temp/2+sqrt((Sigma_M_temp/2)^2+tau^2))/10^6;
        Sigma_min_2 = (Sigma_M_temp/2-sqrt((Sigma_M_temp/2)^2+tau^2))/10^6;


    case 'C'
        SecArea = (2*X*Ty+(Y-2*Ty)*Tx)/100*10^(-4);
        Ixx = ((Tx*(Y-2*Ty)^3/12)+2*((Ty*Ty*Ty*X/12)+(Ty*X*(Y-Ty)*(Y-Ty)/4)))/10000*10^(-8);
%         ElstMod = 20*Ixx/Y;
        Sigma_T = F_a/SecArea;
        Sigma_M = (F_p*L)*(X/2)/Ixx*10^(-3);
        Sigma_max = (Sigma_T+Sigma_M)/10^6;
        Sigma_min = 0;

        Q = (Ty*Y)*(X/2-Ty/2)*10^(-9);
        tau = F_p*Q/Ixx/Tx*10^(3);
        Sigma_M_temp = (F_p*L)*(X/2-Ty)/(210*10^(9)*Ixx)*10^(-3)+Sigma_T;
        Sigma_max_2 = (Sigma_M_temp/2+sqrt((Sigma_M_temp/2)^2+tau^2))/10^6;
        Sigma_min_2 = (Sigma_M_temp/2-sqrt((Sigma_M_temp/2)^2+tau^2))/10^6;

    case 'I'
        SecArea = (2*X*Ty+(Y-2*Ty)*Tx)/100*10^(-4);
        Ixx = ((Tx*(Y-2*Ty)^3/12)+2*((Ty*Ty*Ty*X/12)+(Ty*X*(Y-Ty)*(Y-Ty)/4)))*10^(-12);
        
%         ElstMod = 20*Ixx/Y;
        Sigma_T = F_a/SecArea;
        Sigma_M = (F_p*L)*(X/2)/Ixx*10^(-3);
        Sigma_max = (Sigma_T+Sigma_M)/10^6;
        Sigma_min = 0;

        Q = (Ty*Y)*(X/2-Ty/2)*10^(-9);
        tau = F_p*Q/Ixx/Tx*10^(3);
        Sigma_M_temp = (F_p*L)*(X/2-Ty)/Ixx*10^(-3)+Sigma_T;
        Sigma_max_2 = (Sigma_M_temp/2+sqrt((Sigma_M_temp/2)^2+tau^2))/10^6;
        Sigma_min_2 = (Sigma_M_temp/2-sqrt((Sigma_M_temp/2)^2+tau^2))/10^6;

end
end