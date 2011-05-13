%------------------------------------------------------------------------------
% Run this in matlab/octave after creating section.dat files using slice.py
%------------------------------------------------------------------------------
data=load('section.dat');
n = length(data(:,1));

for j=1:n
   x = data(j,1);
   y = data(j,2);
   Cp(j)= data(j,4);

   theta(j) = atan2(y,x);
   Cpexact(j) = 1 - (9/4)*sin(theta(j))^2;
end

theta=abs(theta);
theta=pi - theta;
plot(theta,Cpexact,'-',theta,-data(:,4),'LineWidth',2)
legend('Exact','flo3d')
ylabel('C_p')
xlabel('\theta')
axis([0 pi -1.5 1.5])
