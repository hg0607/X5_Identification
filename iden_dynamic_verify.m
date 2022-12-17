
clear
close all
format short
file = 'IdenData2022_11_29_10_31_29.csv';
g = 9.81;

pos1 = readmatrix(file,'Range','A2:A100000')*360/524288;
vel1 = readmatrix(file,'Range','D2:D100000')*360/524288;
cur1 = readmatrix(file,'Range','G2:G100000');

pos2 = readmatrix(file,'Range','B2:B100000')*360/524288;
vel2 = readmatrix(file,'Range','E2:E100000')*360/524288;
cur2 = readmatrix(file,'Range','H2:H100000');

pos3 = readmatrix(file,'Range','C2:C100000')*360/524288;
vel3 = readmatrix(file,'Range','F2:F100000')*360/524288;
cur3 = readmatrix(file,'Range','I2:I100000');
n = 3;
m = length(pos1);
t=0.001*(1:1:n);

fc = 5;
fs = 1000;
[b,a] = butter(2,fc/(fs/2));
x3 = filtfilt(b,a,vel3);
y3 = filtfilt(b,a,cur3);
x2 = filtfilt(b,a,vel2);
y2 = filtfilt(b,a,cur2);
x1 = filtfilt(b,a,vel1);
y1 = filtfilt(b,a,cur1);
acc1 = gradient(x1)*1000/1.18;
acc2 = gradient(x2)*1000/1.18;
acc3 = gradient(x3)*1000/1.18;

% Pb_static = [fc_1; fo_1; ...
%             (l_2x + m_3/4)*g; (l_2y+l_3z)*g; fc_2; fo_2; ...
%             l_3x*g; l_3y*g; fc_3; fo_3];
% Pb_dynamic = [Ia_1; fv_1; fc_1; fo_1;  ...
%             (l_2x + m_3/4)*g; (l_2y+l_3z)*g; Ia_2; fv_2; fc_2; fo_2;  ...
%             l_3x*g; l_3y*g; Ia_3; fv_3; fc_3; fo_3 ];
p_static = 10;
Hb_static = zeros(n*m,p_static);
p_dynamic = 16;
Hb_dynamic = zeros(n*m,p_dynamic);
y = zeros(n*m,1);
q1 = deg2rad(pos1); q2 = deg2rad(pos2+90); q3 = deg2rad(pos3);
dq1 = x1; dq2 = x2; dq3 = x3;
for i=1:m
    Hb_static(3*i-2,:) = [tanh(5*dq1(i)), 1,...
           0, 0, 0, 0,...
           0, 0, 0, 0];
       
	Hb_static(3*i-1,:) = [0, 0,...
           -cos(q2(i)), sin(q2(i)), tanh(5*dq2(i)),1,...
           -cos(q2(i))*cos(q3(i)), cos(q2(i))*sin(q3(i)), 0, 0 ];

	Hb_static(3*i,:) = [0, 0,...
           0, 0, 0, 0,...
           sin(q2(i))*sin(q3(i)), sin(q2(i))*cos(q3(i)), tanh(5*dq3(i)),1 ];

    Hb_dynamic(3*i-2,:) = [acc1(i), dq1(i), tanh(5*dq1(i)), 1,...
           0, 0, 0, 0, 0, 0, ...
           0, 0, 0, 0, 0, 0];
       
	Hb_dynamic(3*i-1,:) = [0, 0, 0, 0, ...
           -cos(q2(i)), sin(q2(i)), acc2(i), dq2(i), tanh(5*dq2(i)),1,...
           -cos(q2(i))*cos(q3(i)), cos(q2(i))*sin(q3(i)), 0, 0, 0, 0 ];

	Hb_dynamic(3*i,:) = [0, 0, 0, 0, ...
           0, 0, 0, 0, 0, 0, ...
           sin(q2(i))*sin(q3(i)), sin(q2(i))*cos(q3(i)), acc3(i), dq3(i), tanh(5*dq3(i)),1 ];
       
    y(3*i-2) = y1(i);
    y(3*i-1) = y2(i);
    y(3*i) = y3(i);
end
load x_ls_static.mat 
y_ls_static = Hb_static*x_ls_static;
cond_Hb_static = cond(Hb_static)
figure('name','静力学拟合')
subplot(3,1,1)
plot(y1,'k')
hold on
plot(y_ls_static(1:3:end),'r')
legend('实际值','拟合值')
mse1 = mse(y1 - y_ls_static(1:3:end))
subplot(3,1,2)
plot(y2,'k')
hold on
plot(y_ls_static(2:3:end),'r')
legend('实际值','拟合值')
mse2 = mse(y2 - y_ls_static(2:3:end))
subplot(3,1,3)
plot(y3,'k')
hold on
plot(y_ls_static(3:3:end),'r')
legend('实际值','拟合值')
mse3 = mse(y3 - y_ls_static(3:3:end))

cond_Hb_dynamic = cond(Hb_dynamic)
load x_ls_dynamic.mat 
% x_ls_var = diag(inv(Hb_dynamic'*Hb_dynamic))
y_ls_dynamic = Hb_dynamic*x_ls_dynamic;
figure('name','动力学拟合')
subplot(3,1,1)
plot(y1,'k')
hold on
plot(y_ls_dynamic(1:3:end),'r')
legend('实际值','拟合值')
mse1 = mse(y1 - y_ls_dynamic(1:3:end))
subplot(3,1,2)
plot(y2,'k')
hold on
plot(y_ls_dynamic(2:3:end),'r')
legend('实际值','拟合值')
mse2 = mse(y2 - y_ls_dynamic(2:3:end))
subplot(3,1,3)
plot(y3,'k')
hold on
plot(y_ls_dynamic(3:3:end),'r')
legend('实际值','拟合值')
mse3 = mse(y3 - y_ls_dynamic(3:3:end))

% tau1 = fc1*sign(vel1) + fv1*vel1 + fo1 + Ia1*acc1 % pos单位：度 vel单位：度/秒
% tau2 = -(l_2x + m_3/4)*g*cospos2 + (l_2y + l_3z)*g*sinpos2 -
% l_3x*g*cospos2.*cospos3 + l_3y*g*cospos2.*sinpos3 + fc2*sign(vel2) + fv2*vel2
% + fo2 + Ia2*acc2
% tau3 = l_3x*g*sinpos2.*sinpos3 + l_3y*g*sinpos2.*cospos3 + fc3*sign(vel3) +
% fv3*vel3 + fo3 + Ia3*acc3
