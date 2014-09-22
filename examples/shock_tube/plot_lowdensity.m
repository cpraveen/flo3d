data = load('line.dat');
exact= load('lowdensity.dat');

figure(1)
plot(data(:,1),data(:,2),'o',exact(:,1),exact(:,2),'-')
legend('flo3d','Exact')
title('Density')
axis([0 1 0 1.2])

figure(2)
plot(data(:,1),data(:,3),'o',exact(:,1),exact(:,3),'-')
legend('flo3d','Exact')
title('Velocity')
axis([0 1 -2.5 2.5])

figure(3)
plot(data(:,1),data(:,4),'o',exact(:,1),exact(:,4),'-')
legend('flo3d','Exact')
title('Pressure')
axis([0 1 0 0.5])
