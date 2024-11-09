function [t, y] = rk4(f, a, b, ya, nstep)

% Get the dimension of the system of odes
m = max(size(ya));
% Initialize the t & y lists
t = zeros(1,nstep+1);
y = zeros(m,nstep+1);
% Get the time step (or h in notes)
dt = (b - a)/nstep;
% Set the initial values of the system
t(1) = a;
y(:,1) = ya;

%  0  |  0   0   0   0
%  1  | 1/2  0   0   0
% 1/2 |  0  1/2  0   0
% 1/2 |  0   0  1/2  0
% ----+----------------- 
%     | 1/6 1/3 1/3 1/6
% Use standard RK4 Scheme
for j = 1:nstep
    % Update t
    t(j + 1) = t(j) + dt;
    % RK2 Scheme
    k1 = f(t(j), y(:,j));
    k2 = f(t(j) + dt/2, y(:,j) + 0.5*dt*k1);
    k3 = f(t(j) + dt/2, y(:,j) + 0.5*dt*k2);
    k4 = f(t(j) + dt, y(:,j) + dt*k3);
    y(:,j+1) = y(:,j) + dt*(k1/6 + k2/3 + k3/3 + k4/6);
end

end