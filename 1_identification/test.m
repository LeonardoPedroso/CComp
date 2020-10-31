
clear;
tic
Ts = 0.01; % (s)
t = [0:Ts:100]';
u = sin((2*pi/5)*t);
y = simOpenLoop(t,u,Ts);
toc
figure;
hold on;
plot(t,u);
plot(t,y);
hold off;

