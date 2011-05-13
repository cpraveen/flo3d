%------------------------------------------------------------------------------
% Run this in matlab/octave after creating section?.dat files using slice.py
%------------------------------------------------------------------------------
figure(1)
data=load('section0.dat');
e=load('exp/cp90b.ex');

c = 1.928 - 0.76;
x = (data(:,1) - min(data(:,1))) / c;
plot(e(:,1),-e(:,2),'o',x,data(:,4),'-','LineWidth',2)
axis([0,1,-0.8, 0.6])
xlabel('x/L')
ylabel('-C_p')
title('90 deg.')
%------------------------------------------------------------------------------
figure(2)
data=load('section1.dat');
e=load('exp/cp15b.ex');

c = 1.928 - 0.76;
x = (data(:,1) - min(data(:,1))) / c;
plot(e(:,1),-e(:,2),'o',x,data(:,4),'-','LineWidth',2)
axis([0,1,-0.8, 0.6])
xlabel('x/L')
ylabel('-C_p')
title('15 deg.')
%------------------------------------------------------------------------------
figure(3)
data=load('section2.dat');
e=load('exp/cp25.ex');

c = max(data(:,1)) - min(data(:,1));
x = (data(:,1) - min(data(:,1))) / c;
plot(e(:,1),-e(:,2),'o',x,data(:,4),'-','LineWidth',2)
xlabel('x/c')
ylabel('-C_p')
title('25% span')
%------------------------------------------------------------------------------
figure(4)
data=load('section3.dat');
e=load('exp/cp40.ex');

c = max(data(:,1)) - min(data(:,1));
x = (data(:,1) - min(data(:,1))) / c;
plot(e(:,1),-e(:,2),'o',x,data(:,4),'-','LineWidth',2)
xlabel('x/c')
ylabel('-C_p')
title('40% span')
%------------------------------------------------------------------------------
figure(5)
data=load('section4.dat');
e=load('exp/cp60.ex');

c = max(data(:,1)) - min(data(:,1));
x = (data(:,1) - min(data(:,1))) / c;
plot(e(:,1),-e(:,2),'o',x,data(:,4),'-','LineWidth',2)
xlabel('x/c')
ylabel('-C_p')
title('60% span')
%------------------------------------------------------------------------------
figure(6)
data=load('section5.dat');
e=load('exp/cp75.ex');

c = max(data(:,1)) - min(data(:,1));
x = (data(:,1) - min(data(:,1))) / c;
plot(e(:,1),-e(:,2),'o',x,data(:,4),'-','LineWidth',2)
xlabel('x/c')
ylabel('-C_p')
title('75% span')
%------------------------------------------------------------------------------
figure(7)
data=load('section6.dat');
e=load('exp/cp85.ex');

c = max(data(:,1)) - min(data(:,1));
x = (data(:,1) - min(data(:,1))) / c;
plot(e(:,1),-e(:,2),'o',x,data(:,4),'-','LineWidth',2)
xlabel('x/c')
ylabel('-C_p')
title('85% span')
%------------------------------------------------------------------------------
figure(8)
data=load('section7.dat');
e=load('exp/cp925.ex');

c = max(data(:,1)) - min(data(:,1));
x = (data(:,1) - min(data(:,1))) / c;
plot(e(:,1),-e(:,2),'o',x,data(:,4),'-','LineWidth',2)
xlabel('x/c')
ylabel('-C_p')
title('92.5% span')
