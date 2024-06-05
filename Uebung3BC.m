clear all
close all
clc
addpath("./Funktion")

w_E = [0;0;7.292115e-5];
Omega_iee = [0 -w_E(3) 0
             w_E(3) 0 0
             0 0 0];

% 读取IMU数据 ax ay az(m/s^2) gx gy gz(rad/s) 0.02s
IMU.raw = csvread('INav-Uebung03-A1-IMU.csv', 1, 0);
t = IMU.raw(:,1);
IMUdata(:,1:3) = IMU.raw(:,2:4);
IMUdata(:,4:6) = IMU.raw(:,5:7);
REF.raw = csvread('INav-Uebung03-A1-REF.csv', 1, 0);%0.2s
refp = REF.raw(:,2:4);
refv = REF.raw(:,5:7)';
t_ref = REF.raw(:,1);
%角度转弧度 deg2rad
refRPY = deg2rad(REF.raw(:,8:10));
refp(:,1:2) = deg2rad(refp(:,1:2)); 

%Ref Position in e-system
wgs84 = wgs84Ellipsoid('meter');
[xRef(:),yRef(:),zRef(:)] = geodetic2ecef(refp(:,1),refp(:,2),refp(:,3),wgs84);

%DCM t0 berechnen
Cne_0 = C(3,-refp(1,2))*C(2,refp(1,1)+pi/2);
Cbn_0 = C(3,-refRPY(1,3))*C(2,-refRPY(1,2))*C(1,-refRPY(1,1));
Cpe_0 = Cne_0*Cbn_0;
Omega_iep = inv(Cpe_0)*Omega_iee*Cpe_0;
w_ieP(1,1) = Omega_iep(3,2);
w_ieP(2,1) = Omega_iep(1,3);
w_ieP(3,1) = Omega_iep(2,1);
refv_e(:,1) = Cne_0*refv(:,1);           % v in t0 in e ; refv in n                                     

%DCM t1 berechnen
Cne_1 = C(3,-refp(2,2))*C(2,refp(2,1)+pi/2);
Cbn_1 = C(3,-refRPY(2,3))*C(2,-refRPY(2,2))*C(1,-refRPY(2,1));
Cpe_1 = Cne_1*Cbn_1;
Omega_iep = inv(Cpe_1)*Omega_iee*Cpe_1;
w_ieP(1,1) = Omega_iep(3,2);
w_ieP(2,1) = Omega_iep(1,3);
w_ieP(3,1) = Omega_iep(2,1);
refv_e(:,2) = Cne_1*refv(:,2);           % v in t1 in e ; refv in n                                       

% ref pos 
x(:,1) = xRef(:);
x(:,2) = yRef(:);
x(:,3) = zRef(:);
                                                

%Quaternion t0 im n-system berechnen 
qt0 = compact(quaternion((Cbn_0),'rotmat','frame'));

%Quaternion t1 im n-system berechnen
qt1 = compact(quaternion((Cbn_1),'rotmat','frame'));
q_n(1,:) = qt0;
q_n(2,:) = qt1;

%Quaternion t0 im e-system berechnen 
qt0 = compact(quaternion((Cpe_0),'rotmat','frame'));

%Quaternion t1 im e-system berechnen
qt1 = compact(quaternion((Cpe_1),'rotmat','frame'));
q_e(1,:) = qt0;
q_e(2,:) = qt1;

%% Aufgabe 1.1 n system
Rk3_temp = [refp(1,:),refv(:,1)',q_n(1,:)];%deg NED q
for i = 2:length(t)
    delta_t = abs(t(i) - t(i-1));
    if i == 2
        x0 = [IMUdata(i-1, 1:6); IMUdata(i, 1:6)];
    else
        x0 = [IMUdata(i-1:i, 1:3) IMUdata(i-1:i, 4:6)];
    end
    Rk3_temp(i,:) = Heun(@TimeDerivativePosVelAtt_n, Rk3_temp(i-1,:)', x0', delta_t);

      
end
RK3_Pos = Rk3_temp(:,1:3);
RK3_Vel = Rk3_temp(:,4:6);
RK3_q = Rk3_temp(:,7:10);

RK3_Pos_e = lla2ecef([rad2deg(Rk3_temp(:,1:2)) Rk3_temp(:, 3)]);

figure
plot3(x(:,1),x(:,2),x(:,3))
hold on
plot3(RK3_Pos_e(:,1),RK3_Pos_e(:,2),RK3_Pos_e(:,3))
xlabel("[m]");ylabel("[m]");zlabel("[m]");
title('Aufgabe 1 3D Plot for all');
legend('REF', 'V0304', 'RK3');

figure
plot(x(:,2),x(:,3));
hold on 
plot(RK3_Pos_e(:,2),RK3_Pos_e(:,3));
title('Aufgabe 1 in n-system Plot for YZ');
xlabel("[m]");ylabel("[m]");
legend('REF', 'V0304', 'RK3');

figure
plot(t_ref,x(:,3));
hold on
plot(t,RK3_Pos_e(:,3));


%% Aufgabe 1.2 e-system

Rk3_temp = [x(1,:),refv_e(:,1)',q_e(1,:)];
for i = 2:length(t)
    delta_t = abs(t(i) - t(i-1));
    if i == 2
        x0 = [IMUdata(i-1, 1:6); IMUdata(i, 1:6)];
    else
        x0 = [IMUdata(i-1:i, 1:3) IMUdata(i-1:i, 4:6)];
    end
    Rk3_temp(i,:) = Heun(@TimeDerivativePosVelAtt_e, Rk3_temp(i-1,:)', x0', delta_t);

         
end
RK3_Pos = Rk3_temp(:,1:3);
RK3_Vel = Rk3_temp(:,4:6);
RK3_q = Rk3_temp(:,7:10);

for i = 1:length(RK3_q)
    RK3_rotm{i,1} = rotmat(quaternion(RK3_q(i,:)),'frame');
end

figure
plot3(x(:,1),x(:,2),x(:,3))
hold on
plot3(RK3_Pos(:,1),RK3_Pos(:,2),RK3_Pos(:,3))
xlabel("[m]");ylabel("[m]");zlabel("[m]");
title('Aufgabe 1.2 3D Plot for all');
legend('REF','RK3');

figure
plot(x(:,2),x(:,3));
hold on 

plot(RK3_Pos(:,2),RK3_Pos(:,3));
title(['Aufgabe 1.2 YZ-ECEF Plot']);
xlabel("[m]");ylabel("[m]");
legend('REF', 'RK3');

figure
plot(t_ref,x(:,3));
hold on
plot(t,RK3_Pos(:,3));

%两种方法对比
figure
plot(x(:,2),x(:,3));
hold on 
plot(RK3_Pos(:,2),RK3_Pos(:,3));
plot(RK3_Pos_e(:,2),RK3_Pos_e(:,3));
title('Diff zwischen e- und n-system Plot for YZ');
xlabel("[m]");ylabel("[m]");
legend('REF', 'e', 'n');

figure
plot(t_ref,x(:,3));
hold on
plot(t,RK3_Pos(:,3));
plot(t,RK3_Pos_e(:,3));

title('Diff zwischen e- und n-system Plot for h');
xlabel("[m]");ylabel("[m]");
legend('REF', 'e', 'n');

diffz = RK3_Pos_e(:,3)-RK3_Pos(:,3);
figure
subplot(3,1,1)
plot(diffz);
title('Diff zwischen e- und n-system Plot for h');

diffx = RK3_Pos_e(:,1)-RK3_Pos(:,1);
subplot(3,1,2)
plot(diffx);
title('Diff zwischen e- und n-system Plot for x');

diffy = RK3_Pos_e(:,2)-RK3_Pos(:,2);
subplot(3,1,3)
plot(diffy);
title('Diff zwischen e- und n-system Plot for y');
%% Aufgabe 2
IMU = readtable("./INav-Uebung03-A2-IMU.csv",'VariableNamingRule','preserve');
REF = readtable("./INav-Uebung03-A2-REF.csv",'VariableNamingRule','preserve');
GNSS = readtable("./INav-Uebung03-A2-GNSS.csv",'VariableNamingRule','preserve');

Mes.t = table2array(IMU(:,1));
Mes.acc = table2array(IMU(:,2:4));
Mes.omega = table2array(IMU(:,5:7));

Ref.t = table2array(REF(:,1));
Ref.Pos = table2array(REF(:,2:4));
Ref.Vel = table2array(REF(:,5:7));

gnss.t = table2array(GNSS(:,1));
gnss.Pos = table2array(GNSS(:,2:4));
gnss.Vel = table2array(GNSS(:,5:7));
%% Standardabweichung berechnen %% 这里是有问题的
acc_std = std(Mes.acc);
omega_std = std(Mes.omega);
pos_std = std(gnss.Pos);
vel_std = std(gnss.Vel);
pos_std_R = std(Ref.Pos);
vel_std_R = std(Ref.Vel);
%% Analysieren
figure
plot3(gnss.Pos(:,1),gnss.Pos(:,2),gnss.Pos(:,3))
hold on
plot3(Ref.Pos(:,1),Ref.Pos(:,2),Ref.Pos(:,3))
figure
plot3(gnss.Vel(:,1),gnss.Vel(:,2),gnss.Vel(:,3))
hold on
plot3(Ref.Vel(:,1),Ref.Vel(:,2),Ref.Vel(:,3))

figure
plot3(Ref.Pos(:,1)-gnss.Pos(:,1),Ref.Pos(:,2)-gnss.Pos(:,2),Ref.Pos(:,3)-gnss.Pos(:,3))

titles = ["Acc-X","Acc-Y","Acc-Z","Gyro-X","Gyro-Y","Gyro-Z"];
ylabels = ["m/s^2","m/s^2","m/s^2","rad/s","rad/s","rad/s"];
figure,
for i = 1:3
    subplot(3,1,i)
    plot(Mes.t,Mes.acc(:,i))
    xlabel('t [s]'),ylabel(ylabels(i)),title(titles(i))
    Std.acc(i,1) = std(Mes.acc(:,i));
end

figure,
for i = 1:3
    subplot(3,1,i)
    plot(Mes.t,Mes.omega(:,i))
    xlabel('t [s]'),ylabel(ylabels(i+3)),title(titles(i+3))
    Std.omega(i,1) = std(Mes.omega(:,i));
end

titles = ["Pos-X","Pos-Y","Pos-Z","Vel-X","Vel-Y","Vel-Z"];
ylabels = ["m","m","m","m/s","m/s","m/s"];
figure,
for i = 1:3
    subplot(3,1,i)
    plot(Ref.t,Ref.Pos(:,i))
    hold on
    plot(gnss.t,gnss.Pos(:,i))
    legend('REF','GNSS')
    xlabel('t [s]'),ylabel(ylabels(i)),title(titles(i))
    Std.Pos(i,1) = std(Ref.Pos(:,i));
end
figure,
for i = 1:3
    subplot(3,1,i)
    plot(Ref.t,Ref.Vel(:,i))
    hold on
    plot(gnss.t,gnss.Vel(:,i))
    legend('REF','GNSS')
    xlabel('t [s]'),ylabel(ylabels(i+3)),title(titles(i+3))
    Std.Vel(i,1) = std(Ref.Vel(:,i));
end





