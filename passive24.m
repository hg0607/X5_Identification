clear
close all
format short
file = 'X5Data2022_12_12_18_57_40.csv';
syms L_1xx L_1xy L_1xz L_1yy L_1yz L_1zz l_1x l_1y l_1z m_1 Ia_1 fv_1 fc_1 fo_1;
syms L_2xx L_2xy L_2xz L_2yy L_2yz L_2zz l_2x l_2y l_2z m_2 Ia_2 fv_2 fc_2 fo_2;
syms L_3xx L_3xy L_3xz L_3yy L_3yz L_3zz l_3x l_3y l_3z m_3 Ia_3 fv_3 fc_3 fo_3;

g = 9.81;
time = readmatrix(file,'Range','E2:E100000');
n = length(time);
msec = zeros(n,1);
for i = 1:1:n-1
    if(time(i)>time(i+1))
        msec(i) = time(i+1)+1000-time(i);
    else
        msec(i) = time(i+1)-time(i);
    end    
end
mean(msec)
std(msec)

pos1 = readmatrix(file,'Range','F16000:F100000');
pos2 = readmatrix(file,'Range','G16000:G100000');
pos3 = readmatrix(file,'Range','H16000:H100000');

vel1 = readmatrix(file,'Range','I16000:I100000');
vel2 = readmatrix(file,'Range','J16000:J100000');
vel3 = readmatrix(file,'Range','K16000:K100000');

cur1 = readmatrix(file,'Range','L16000:L100000');
cur2 = readmatrix(file,'Range','M16000:M100000');
cur3 = readmatrix(file,'Range','N16000:N100000');

n = 3;
m = length(pos1);
fc = 5;
fs = 500;
[b,a] = butter(2,fc/(fs/2));
x3 = filtfilt(b,a,vel3);
y3 = filtfilt(b,a,cur3);
x2 = filtfilt(b,a,vel2);
y2 = filtfilt(b,a,cur2);
x1 = filtfilt(b,a,vel1);
y1 = filtfilt(b,a,cur1);
acc1 = gradient(x1)*1000/5;
acc2 = gradient(x2)*1000/5;
acc3 = gradient(x3)*1000/5;

% Pb_dynamic = [Ia_1; fv_1; fc_1; fo_1;  ...
%             (l_2x + m_3/4)*g; (l_2y+l_3z)*g; Ia_2; fv_2; fc_2; fo_2;  ...
%             l_3x*g; l_3y*g; Ia_3; fv_3; fc_3; fo_3 ];
r = 100;
pb = 24;
Hb_dynamic = zeros(n*m,pb);
H = zeros(m,n*pb);
y = [y1;y2;y3];
q1 = deg2rad(pos1); q2 = deg2rad(pos2+90); q3 = deg2rad(pos3);
dq1 = deg2rad(x1); dq2 = deg2rad(x2); dq3 = deg2rad(x3);
q = [q1 q2 q3];
dq = [dq1 dq2 dq3];
ddq = [acc1 acc2 acc3];

data = importdata('E:\src\SymPyBotics\H_mat.txt');
eval(cell2mat(data(1)));
len = length(data);
for i=2:len
    eval(cell2mat(data(i)));
end

for j=1:pb
    Hb_dynamic(:,j) = [H(:,3*j-2); H(:,3*j-1); H(:,3*j)];
    
end
cond(Hb_dynamic)

load E:\src\SymPyBotics\x_ls_dynamic.mat 
% x_ls_var = diag(inv(Hb_dynamic'*Hb_dynamic))
y_ls_dynamic = Hb_dynamic*x_ls_dynamic;
figure('name','动力学拟合')
subplot(3,1,1)
plot(y1,'k')
hold on
plot(y_ls_dynamic(1:m),'r')
legend('实际值','拟合值')
mse1 = mse(y1 - y_ls_dynamic(1:m))
subplot(3,1,2)
plot(y2,'k')
hold on
plot(y_ls_dynamic(m+1:2*m),'r')
legend('实际值','拟合值')
mse2 = mse(y2 - y_ls_dynamic(m+1:2*m))
subplot(3,1,3)
plot(y3,'k')
hold on
plot(y_ls_dynamic(2*m+1:3*m),'r')
legend('实际值','拟合值')
mse3 = mse(y3 - y_ls_dynamic(2*m+1:3*m))

% tau1 = fc1*sign(vel1) + fv1*vel1 + fo1 + Ia1*acc1 % pos单位：度 vel单位：度/秒
% tau2 = -(l_2x + m_3/4)*g*cospos2 + (l_2y + l_3z)*g*sinpos2 -
% l_3x*g*cospos2.*cospos3 + l_3y*g*cospos2.*sinpos3 + fc2*sign(vel2) + fv2*vel2
% + fo2 + Ia2*acc2
% tau3 = l_3x*g*sinpos2.*sinpos3 + l_3y*g*sinpos2.*cospos3 + fc3*sign(vel3) +
% fv3*vel3 + fo3 + Ia3*acc3

function y = Matrix(x)
    y = x;
end