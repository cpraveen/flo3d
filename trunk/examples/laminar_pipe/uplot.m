data = load('u.dat');

plot(data(:,1),data(:,2),'-',data(:,1),data(:,3),'--','LineWidth',1.5)
legend('flo3d','Exact')
xlabel('y')
ylabel('u')
title('Axial velocity')
