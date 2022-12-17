
clear
close all
format short
file = 'trajectory_record(1).csv';
g = 9.81;

pos1 = readmatrix(file,'Range','A1:A100000')*360/524288;
vel1 = readmatrix(file,'Range','B1:B100000')*360/524288;
cur1 = readmatrix(file,'Range','C1:C100000');

pos2 = readmatrix(file,'Range','D1:D100000')*360/524288;
vel2 = readmatrix(file,'Range','E1:E100000')*360/524288;
cur2 = readmatrix(file,'Range','F1:F100000');

pos3 = readmatrix(file,'Range','G1:G100000')*360/524288;
vel3 = readmatrix(file,'Range','H1:H100000')*360/524288;
cur3 = readmatrix(file,'Range','I1:I100000');
n = 3;
m = length(pos1);
t=0.001*(1:1:m);

fc = 5;
fs = 1000;
[b,a] = butter(2,fc/(fs/2));
x3 = filtfilt(b,a,vel3);
y3 = filtfilt(b,a,cur3);
x2 = filtfilt(b,a,vel2);
y2 = filtfilt(b,a,cur2);
x1 = filtfilt(b,a,vel1);
y1 = filtfilt(b,a,cur1);

% Pb = [fc_1; fo_1;...
%       (l_2x + m_3/4)*g; (l_2y+l_3z)*g; fc_2; fo_2;...
%       l_3x*g; l_3y*g; fc_3; fo_3];
p = 10;
Hb = zeros(n*m,p);
y = zeros(n*m,1);
q1 = deg2rad(pos1); q2 = deg2rad(pos2+90); q3 = deg2rad(pos3);
dq1 = x1; dq2 = x2; dq3 = x3;
for i=1:m
    Hb(3*i-2,:) = [tanh(5*dq1(i)), 1,...
           0, 0, 0, 0,...
           0, 0, 0, 0];
       
	Hb(3*i-1,:) = [0, 0,...
           -cos(q2(i)), sin(q2(i)), tanh(5*dq2(i)), 1,...
           -cos(q2(i))*cos(q3(i)), cos(q2(i))*sin(q3(i)), 0, 0 ];

	Hb(3*i,:) = [0, 0,...
           0, 0, 0, 0,...
           sin(q2(i))*sin(q3(i)), sin(q2(i))*cos(q3(i)), tanh(5*dq3(i)),1 ];

       
    y(3*i-2) = y1(i);
    y(3*i-1) = y2(i);
    y(3*i) = y3(i);
end
cond_Hb = cond(Hb)
x_ls_static =(Hb'*Hb)\Hb'*y
save 'x_ls_static.mat' x_ls_static
x_ls_var = diag(inv(Hb'*Hb))
y_ls = Hb*x_ls_static;
figure

subplot(3,1,1)
plot(y1,'k')
hold on
plot(y_ls(1:3:end),'r')
legend('实际值','拟合值')
mse1 = mse(y1 - y_ls(1:3:end))

subplot(3,1,2)
plot(y2,'k')
hold on
plot(y_ls(2:3:end),'r')
legend('实际值','拟合值')
mse2 = mse(y2 - y_ls(2:3:end))

subplot(3,1,3)
plot(y3,'k')
hold on
plot(y_ls(3:3:end),'r')
legend('实际值','拟合值')
mse3 = mse(y3 - y_ls(3:3:end))

% tau1 = fc1*tanh(vel1) + fo1 % pos单位：度 vel单位：度/秒 
% tau2 = -(l_2x + m_3/4)*g*cospos2 + (l_2y + l_3z)*g*sinpos2 - l_3x*g*cospos2.*cospos3 + l_3y*g*cospos2.*sinpos3 + fc2*tanh(vel2) + fo2 
% tau3 = l_3x*g*sinpos2.*sinpos3 + l_3y*g*sinpos2.*cospos3 + fc3*tanh(vel3) + fo3
