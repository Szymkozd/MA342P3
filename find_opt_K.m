function u = find_opt_K(K,t_vec_end)

% Parameters, BCs, ICs
H_rho = 0.1;
L = 2;
dx = 0.01;
dt = 0.028;

x_vec = 0:dx:L;
t_vec = 0:dt:t_vec_end;
u = zeros(length(x_vec),length(t_vec));
u(:,2) = (exp(-10*(x_vec-(L/2)).^2)-exp(-10*(L/2)^2))*dt;

% Future Displacement
i = 2:1:(length(x_vec)-1);
for j = 2:1:length(t_vec)
    u(i,j+1) = ((H_rho/dx^2)*(u(i+1,j)+u(i-1,j)-2*u(i,j)) - (1/dt^2)*u(i,j-1) - ((-2/dt^2)-(K/dt))*u(i,j))/((K/dt)+(1/dt^2));
end

end