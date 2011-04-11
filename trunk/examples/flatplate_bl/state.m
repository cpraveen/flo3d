format compact
mach_inf = 0.25
rho_inf = 0.1
T_ref = 273.15
mu_ref = 1.716e-5
T_inf = 300
mu_inf = mu_ref * (T_inf/T_ref)^1.5 * (T_ref + 110.4)/(T_inf + 110.4)
gamma=1.4
L = 0.3

R = 287.0;
p_inf = rho_inf * R * T_inf
c_inf = sqrt(gamma * p_inf / rho_inf)
u_inf = c_inf * mach_inf
Re = rho_inf * u_inf * L / mu_inf

% BL thickness at end of plate
d = 4.5*L/sqrt(Re)
eta=50;
y=sqrt(2*L*L/Re)*eta
