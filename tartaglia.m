clear all;
close all;
clc;

tht = 30;
v = 20;
g = 9.81;
m = 6e-3;
d = 70e-3;

% Calculate terminal velocity
vt = sqrt(2 * m * g / (1.225 * 0.1 * pi * d^2 * 0.25));

% Initial velocities
vx0 = v * cosd(tht);
vy0 = v * sind(tht);

% Time vector
ts = 0:0.001:10;

% Differential equations
fun1 = @(t, u) [-g * (sqrt(u(1)^2 + u(2)^2) * u(1) / vt^2); -g * (1 + (sqrt(u(1)^2 + u(2)^2) * u(2) / vt^2))];
fun2 = @(t, u) [-g * (sqrt(u(1)^2 + u(2)^2) * u(1) / vt^2); -g * (1 - (sqrt(u(1)^2 + u(2)^2) * u(2) / vt^2))];

% Solve differential equations
[t, s] = ode45(fun1, ts, [vx0, vy0]);
[val, idx] = min(s(:, 2) >= 0);
[t1, s1] = ode45(fun2, ts(1, idx-1:end), [s(idx-1, 1), s(idx-1, 2)]);

% Extract ascending and descending phases
u_as = s(1:idx-1, 1);
v_as = s(1:idx-1, 2);
u_des = s1(:, 1);
v_des = s1(:, 2);

% Calculate total horizontal and vertical distances
dx_as = trapz(ts(1, 1:idx-1), s(1:idx-1, 1));
dy_as = trapz(ts(1, 1:idx-1), s(1:idx-1, 2));

% Find the point where the projectile returns to the ground
for i = 2:size(s1, 1)
    ry = trapz(ts(1, 1:i), s1(1:i, 2));
    rx = trapz(ts(1, 1:i), s1(1:i, 1));
    if ry + dy_as <= 0
        break;
    end
end

% Calculate x and y positions for both phases
x_as = zeros(idx-1, 1);
y_as = zeros(idx-1, 1);
for k = 2:idx-1
    x_as(k) = trapz(t(1:k, 1), u_as(1:k, 1));
    y_as(k) = trapz(t(1:k, 1), v_as(1:k, 1));
end

x_des = zeros(size(s1, 1), 1);
y_des = zeros(size(s1, 1), 1);
for k = 2:size(t1, 1)
    x_des(k) = trapz(t(1:k, 1), u_des(1:k, 1));
    y_des(k) = trapz(t(1:k, 1), v_des(1:k, 1));
end

% Find the point where the projectile returns to the ground for the descending phase
[vi, ii] = min((y_des + y_as(end)) >= 0);

% Plot the projectile trajectory
plot(x_as, y_as, 'b-', 'LineWidth', 2);
hold on;
plot(x_as(end) + x_des(1:ii), y_as(end) + y_des(1:ii), 'r-', 'LineWidth', 2);
xlabel('Horizontal Distance (m)');
ylabel('Vertical Distance (m)');
title('Projectile Trajectory');
legend('Ascending Phase', 'Descending Phase');
grid on;
axis equal;