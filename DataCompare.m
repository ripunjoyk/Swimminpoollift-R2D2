clc;
Sec_in=input('Section shape\n');
D_in=input('D\n');
B_in=input('B\n');
t_in=input('t\n');
T_in=input('T\n');
Data=Crosssectionproperties;
l=height(Data);
Shortlisted_Data=Data(Data.Section == Sec_in, :);
Shortlisted_Data2=Data(Data.D == D_in & Data.B == B_in & Data.t == t_in & Data.T == T_in,:);
Area_actual=Shortlisted_Data2.A;
Ix_actual=Shortlisted_Data2.Ix;
[Ix_cal,E,Area_cal]=crossectional_analysis(Sec_in,D_in,B_in,t_in,T_in);
if height(Shortlisted_Data2)==0
    disp('The given data did not match with any section in the data table.');
else
    Error_area_p=norm(Area_actual-Area_cal)/Area_actual*100;
    Error_Ix_p=norm(Ix_actual-Ix_cal)/Ix_actual*100;
    disp('The error in area in percentage is');
    disp(Error_area_p);
    disp('The error in Ix in percentage is');
    disp(Error_Ix_p);
end
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
        %         Ixx =  (T*(Y-4*T)^3/6)+ (( ((B-4T)*T^3/3)+(T*(B-4T)*(D-T)^2))/2) + (pi*t^4(405-(3136/pi^2))/108)+(3*pi*T^2)*((9*pi*(Y-4*T)+56*T)/18*pi)^2;
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