clc;
%FS=input('Factor of safety\n');
L=.375;%input('Length in metres\n');
F_a=-5000;%input('Axial force\n');
F_p=160000;%input('Perpendicular force\n');
M=0;%input('Moment (+ for shagging)\n');
Section='I';%input('Section shape\n');
X=206;%input('D\n');
Y=206;%input('B\n');
Tx=7.87;%input('t\n');
Ty=12.6;%input('T\n');
% [Sigma_max,Sigma_min,Sigma_max_In,Sigma_min_In]=Stress_P(Section,X,Y,Tx,Ty,F_a,F_p,M,L);
% if Sigma_max>=max(abs(Sigma_max_In),abs(Sigma_min_In))
%     disp('Max stress occurs at the boundary and is given by :');
%     disp(strcat(num2str(Sigma_max,'%.4f'),' MPa'));
% else
%     disp('Max stress occurs at the point between the web and flange and is given by');
%     disp('Sigma_max');
%     disp(strcat(num2str(Sigma_max_In,'%.4f'),' MPa'));
%     disp('Sigma_min');
%     disp(strcat(num2str(Sigma_min_In,'%.4f'),' MPa'));
% end

[Sigma_max,Sigma_min,Sigma_max_2,Sigma_min_2,Sigma_max_In,Sigma_min_In,Sigma_max_In_2,Sigma_min_In_2]=Stress_P(Section,X,Y,Tx,Ty,F_a,F_p,M,L);


%%

function [Ixx,ElstMod,SecArea] = crossectional_analysis(Section,X,Y,Tx,Ty)
switch Section
    case 'P'
        SecArea = (pi*X*Y-pi*(X-2*Tx)*(Y-2*Ty))/400;
        Ixx = pi*(X*Y*X*Y-(X-2*Tx)^4)/640000;
        ElstMod = 20*Ixx/X;

    case 'B'
        SecArea = 2*Tx*((X-4*Tx)+(Y-4*Tx)+(3*pi*Tx/2))/100;
        Ixx =  ((X*Y*Y*Y)-(X-2*Tx)*(Y-2*Ty)^3)/120000;
        ElstMod = 20*Ixx/Y;

    case 'L'
        SecArea = ((X-Tx)*Ty+(Y-Ty)*Tx+Tx*Ty)/100;
        Ixx = ((Tx*(Y-Ty)^3/12)+(X*Ty^3/12)+(X*Ty*(X-Ty)*(X-Ty)/4))/10000;
        ElstMod = 20*Ixx/Y;

    case 'C'
        SecArea = (2*X*Ty+(Y-2*Ty)*Tx)/100;
        Ixx = ((Tx*(Y-2*Ty)^3/12)+2*((Ty*Ty*Ty*X/12)+(Ty*X*(Y-Ty)*(Y-Ty)/4)))/10000;
        ElstMod = 20*Ixx/Y;

    case 'I'
        SecArea = (2*X*Ty+(Y-2*Ty)*Tx)/100;
        Ixx = (Tx*(Y-2*Ty)^3/12)+2*((Ty*Ty*Ty*X/12)+(Ty*X*(Y-Ty)*(Y-Ty)/4));
        Ixx = Ixx/10000;
        ElstMod = 20*Ixx/Y;
end
end

%% 
function [Sigma_max,Sigma_min,Sigma_max_2,Sigma_min_2,Sigma_max_In,Sigma_min_In,Sigma_max_In_2,Sigma_min_In_2] = Stress_P(Section,X,Y,Tx,Ty,F_a,F_p,M,L)
        Sigma_min = 0;
        Sigma_min_2 = 0;
switch Section
    case 'P'
        [Ixx,ElstMod,SecArea]=crossectional_analysis(Section,X,Y,Tx,Ty);
        SecArea=SecArea*10^-4;
        Ixx=Ixx*10^-8;
        Sigma_T = F_a/SecArea;
        Sigma_M = (F_p*L)*(X/2)/Ixx*10^(-3);
        Sigma_MM = M*(X/2)/Ixx*10^(-3);
        Sigma_max = (Sigma_T+Sigma_M-Sigma_MM)/10^6;
        Sigma_max_2 = (Sigma_T-Sigma_M+Sigma_MM)/10^6;


        Sigma_max_In = Sigma_max;
        Sigma_min_In = 0;
        Sigma_max_In_2 = Sigma_max;
        Sigma_min_In_2 = 0;

    case 'B'
        [Ixx,ElstMod,SecArea]=crossectional_analysis(Section,X,Y,Tx,Ty);
        SecArea=SecArea*10^-4;
        Ixx=Ixx*10^-8;
        Sigma_T = F_a/SecArea;
        Sigma_M = (F_p*L)*(X/2)/Ixx*10^(-3);
        Sigma_MM = M*(X/2)/Ixx*10^(-3);
        Sigma_max = (Sigma_T+Sigma_M-Sigma_MM)/10^6;
        Sigma_max_2 = (Sigma_T-Sigma_M+Sigma_MM)/10^6;
        
        Q = (Ty*Y)*(X/2-Ty/2)*10^(-9);
        tau = F_p*Q/Ixx/Tx/2*10^(3);
        Sigma_MM_2 = M*(X/2-Ty)/Ixx*10^(-3);
        Sigma_M_temp = (F_p*L)*(X/2-Ty)/Ixx*10^(-3)+Sigma_T-Sigma_MM_2;
        Sigma_M_temp_2 = -(F_p*L)*(X/2-Ty)/Ixx*10^(-3)+Sigma_T+Sigma_MM_2;
        Sigma_max_In = (Sigma_M_temp/2+sqrt((Sigma_M_temp/2)^2+tau^2))/10^6;
        Sigma_min_In = (Sigma_M_temp/2-sqrt((Sigma_M_temp/2)^2+tau^2))/10^6;
        Sigma_max_In_2 = (Sigma_M_temp_2/2+sqrt((Sigma_M_temp_2/2)^2+tau^2))/10^6;
        Sigma_min_In_2 = (Sigma_M_temp_2/2-sqrt((Sigma_M_temp_2/2)^2+tau^2))/10^6;
    case 'L'
        [Ixx,ElstMod,SecArea]=crossectional_analysis(Section,X,Y,Tx,Ty);
        SecArea=SecArea*10^-4;
        Ixx=Ixx*10^-8;
        Sigma_T = F_a/SecArea;
        Sigma_M = (F_p*L)*(X/2)/Ixx*10^(-3);
        X=X-((Y*Ty^2/Tx+X^2-Ty^2)/(2*(Ty*Y/Tx-Ty+X)));
        Sigma_MM = M*(X)/Ixx*10^(-3);
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
        [Ixx,ElstMod,SecArea]=crossectional_analysis(Section,X,Y,Tx,Ty);
        SecArea=SecArea*10^-4;
        Ixx=Ixx*10^-8;
        Sigma_T = F_a/SecArea;
        Sigma_M = (F_p*L)*(X/2)/Ixx*10^(-3);
        Sigma_MM = M*(X/2)/Ixx*10^(-3);
        Sigma_max = (Sigma_T+Sigma_M-Sigma_MM)/10^6;
        Sigma_max_2 = (Sigma_T-Sigma_M+Sigma_MM)/10^6;

        Q = (Ty*Y)*(X/2-Ty/2)*10^(-9);
        tau = F_p*Q/Ixx/Tx*10^(3);
        Sigma_MM_2 = M*(X/2-Ty)/Ixx*10^(-3);
        Sigma_M_temp = (F_p*L)*(X/2-Ty)/Ixx*10^(-3)+Sigma_T-Sigma_MM_2;
        Sigma_M_temp_2 = -(F_p*L)*(X/2-Ty)/Ixx*10^(-3)+Sigma_T+Sigma_MM_2;
        Sigma_max_In = (Sigma_M_temp/2+sqrt((Sigma_M_temp/2)^2+tau^2))/10^6;
        Sigma_min_In = (Sigma_M_temp/2-sqrt((Sigma_M_temp/2)^2+tau^2))/10^6;
        Sigma_max_In_2 = (Sigma_M_temp_2/2+sqrt((Sigma_M_temp_2/2)^2+tau^2))/10^6;
        Sigma_min_In_2 = (Sigma_M_temp_2/2-sqrt((Sigma_M_temp_2/2)^2+tau^2))/10^6;

    case 'I'
        [Ixx,ElstMod,SecArea]=crossectional_analysis(Section,X,Y,Tx,Ty);
        SecArea=SecArea*10^-4;
        Ixx=Ixx*10^-8;
        Sigma_T = F_a/SecArea;
        Sigma_M = (F_p*L)*(X/2)/Ixx*10^(-3);
        Sigma_MM = M*(X/2)/Ixx*10^(-3);
        Sigma_max = (Sigma_T+Sigma_M-Sigma_MM)/10^6;
        Sigma_max_2 = (Sigma_T-Sigma_M+Sigma_MM)/10^6;

        Q = (Ty*Y)*(X/2-Ty/2)*10^(-9);
        tau = F_p*Q/Ixx/Tx*10^(3);
        Sigma_MM_2 = M*(X/2-Ty)/Ixx*10^(-3);
        Sigma_M_temp = (F_p*L)*(X/2-Ty)/Ixx*10^(-3)+Sigma_T-Sigma_MM_2;
        Sigma_M_temp_2 = -(F_p*L)*(X/2-Ty)/Ixx*10^(-3)+Sigma_T+Sigma_MM_2;
        Sigma_max_In = (Sigma_M_temp/2+sqrt((Sigma_M_temp/2)^2+tau^2))/10^6;
        Sigma_min_In = (Sigma_M_temp/2-sqrt((Sigma_M_temp/2)^2+tau^2))/10^6;
        Sigma_max_In_2 = (Sigma_M_temp_2/2+sqrt((Sigma_M_temp_2/2)^2+tau^2))/10^6;
        Sigma_min_In_2 = (Sigma_M_temp_2/2-sqrt((Sigma_M_temp_2/2)^2+tau^2))/10^6;

end
A=[Sigma_max,Sigma_min,Sigma_max_2,Sigma_min_2,Sigma_max_In,Sigma_min_In,Sigma_max_In_2,Sigma_min_In_2];
Temporary=max(A);
Temporary2=abs(min(A));
Temporary=max(Temporary,Temporary2);
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