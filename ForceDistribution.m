clc;
clear variables;
%% Mass

mg=250;
Mg=1200;
F_load=mg+Mg;
Fs=1500;

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
temp_angle=acos(Lac/Lce);


%%

theta=linspace(70*pi/180,140*pi/180,100);
phi=(theta-90)*pi/180;
temp=theta-temp_angle;
EK=sqrt(Lck^2+Lce^2-2*Lck*Lce*cos(temp));
psi=asin(EK/(Lce*sin(temp)));

%%
M_load=F_load*Lgj;

Fgx=M_load/Lgh;
FGX=ones(100)*Fgx;
Fgy=F_load*cos(phi);
Fhx=Fgx+F_load*sin(phi);
Fdx=Fhx;
Fcx=Fgx+Fs*cos(temp);
Fcy=Fgy-Fs*sin(temp);

figure(1)
hold on
% plot(FGX);
plot(Fgy);
% plot(Fhx);
% plot(Fdx);
% plot(Fcx);
% plot(Fcy);

