// Clear all variables (helps make script more predictable in behaviour)
clear all
// Unset functions
funcprot(0)
// Integration domain
t0 = 0;
t1 = 100;
// Number of integration points
N  = 100000;
// time vector
t  = linspace(t0, t1, N+1)';

// Here is the Lorenz equation this function represents:
// https://en.wikipedia.org/wiki/Lorenz_system#Overview
function ydot = lorenz(t, y)
    ydot = zeros(3,1);
    P = 10;                               // This is sigma
    r = 28;                               // This is rho
    b = 8/3;                              // This is beta
    ydot(1) = P*(y(2) - y(1));            // This is dx/dt
    ydot(2) = -y(1)*y(3) + r*y(1) - y(2); // This is dy/dt
    ydot(3) = y(1)*y(2) - b*y(3);         // This is dz/dt
endfunction

// Initial conditions
y0 = [1.0; 1.0; 1.0];

// Solution
y  = ode(y0, t0, t', lorenz)';

// Plotting in one figure window
figure(1)
f1=gcf();
f1.background = color('white');

// first subplot showing x vs y
subplot(221)
plot2d(y(:,1),y(:,2),21);
xlabel('x(t)');
ylabel('y(t)');

// second subplot showing x vs z
subplot(222)
plot2d(y(:,1),y(:,3),2);
xlabel('x(t)');
ylabel('z(t)');

// third subplot showing y vs z
subplot(223)
plot2d(y(:,2),y(:,3),15);
xlabel('y(t)');
ylabel('z(t)');

// second figure showing parametric plot of x, y, z
figure(2)
f2 = gcf();
f2.background = color('white');
param3d(y(:,1),y(:,2),y(:,3))
