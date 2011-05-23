qinf = 0.5;
r    = 1.0;
area = 4*pi*r^2;

data=load('force.dat');
iter = data(:,1);

ndata = length(iter);

for j=1:ndata
   cd(j)  = data(j,2) / (qinf * area);
   cl1(j) = data(j,3) / (qinf * area);
   cl2(j) = data(j,4) / (qinf * area);
end

plot(iter, cd, iter, cl)
fprintf(1,'Cd = %e\n', cd(ndata));
fprintf(1,'Cl = %e\n', cl1(ndata));
fprintf(1,'Cl = %e\n', cl2(ndata));
