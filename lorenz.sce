clear all
funcprot(0)
t0 = 0;
t1 = 100;
N  = 100000;
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

subplot(221)
plot2d(y(:,1),y(:,2),21);
xlabel('x(t)');
ylabel('y(t)');

subplot(222)
plot2d(y(:,1),y(:,3),2);
xlabel('x(t)');
ylabel('z(t)');

subplot(223)
plot2d(y(:,2),y(:,3),16);
xlabel('y(t)');
ylabel('z(t)');
