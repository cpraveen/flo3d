%------------------------------------------------------------------------------
% Run this in matlab/octave after creating section?.dat files using slice.py
%------------------------------------------------------------------------------
data=load('section0.dat');
u=load('exp/cp1u.ex');
l=load('exp/cp1l.ex');
e = [u; l];

c = max(data(:,1)) - min(data(:,1));
x = (data(:,1) - min(data(:,1))) / c;
figure(1)
plot(e(:,1),e(:,2),'o',x,data(:,4),'-','LineWidth',2)
xlabel('x/c')
ylabel('-C_p')
title('20% span')
%------------------------------------------------------------------------------
data=load('section1.dat');
u=load('exp/cp2u.ex');
l=load('exp/cp2l.ex');
e = [u; l];

c = max(data(:,1)) - min(data(:,1));
x = (data(:,1) - min(data(:,1))) / c;
figure(2)
plot(e(:,1),e(:,2),'o',x,data(:,4),'-','LineWidth',2)
xlabel('x/c')
ylabel('-C_p')
title('44% span')
%------------------------------------------------------------------------------
data=load('section2.dat');
u=load('exp/cp3u.ex');
l=load('exp/cp3l.ex');
e = [u; l];

c = max(data(:,1)) - min(data(:,1));
x = (data(:,1) - min(data(:,1))) / c;
figure(3)
plot(e(:,1),e(:,2),'o',x,data(:,4),'-','LineWidth',2)
xlabel('x/c')
ylabel('-C_p')
title('65% span')
%------------------------------------------------------------------------------
data=load('section3.dat');
u=load('exp/cp4u.ex');
l=load('exp/cp4l.ex');
e = [u; l];

c = max(data(:,1)) - min(data(:,1));
x = (data(:,1) - min(data(:,1))) / c;
figure(4)
plot(e(:,1),e(:,2),'o',x,data(:,4),'-','LineWidth',2)
xlabel('x/c')
ylabel('-C_p')
title('80% span')
%------------------------------------------------------------------------------
data=load('section4.dat');
u=load('exp/cp5u.ex');
l=load('exp/cp5l.ex');
e = [u; l];

c = max(data(:,1)) - min(data(:,1));
x = (data(:,1) - min(data(:,1))) / c;
figure(5)
plot(e(:,1),e(:,2),'o',x,data(:,4),'-','LineWidth',2)
xlabel('x/c')
ylabel('-C_p')
title('90% span')
%------------------------------------------------------------------------------
data=load('section5.dat');
u=load('exp/cp6u.ex');
l=load('exp/cp6l.ex');
e = [u; l];

c = max(data(:,1)) - min(data(:,1));
x = (data(:,1) - min(data(:,1))) / c;
figure(6)
plot(e(:,1),e(:,2),'o',x,data(:,4),'-','LineWidth',2)
xlabel('x/c')
ylabel('-C_p')
title('95% span')
%------------------------------------------------------------------------------
data=load('section6.dat');
u=load('exp/cp7u.ex');
l=load('exp/cp7l.ex');
e = [u; l];

c = max(data(:,1)) - min(data(:,1));
x = (data(:,1) - min(data(:,1))) / c;
figure(7)
plot(e(:,1),e(:,2),'o',x,data(:,4),'-','LineWidth',2)
xlabel('x/c')
ylabel('-C_p')
title('99% span')
