clear all
a = 0;
b = 100;
N = 100000;
t = linspace(a, b, N+1); 

function ydot = lorenz(t, y)
    ydot = zeros(3,1);
    P = 10;
    r = 28;
    b = 8/3;
    ydot(1) = P*(y(2) - y(1));
    ydot(2) = -y(1)*y(3) + r*y(1) - y(2);
    ydot(3) = y(1)*y(2) - b*y(3);
endfunction

y = ode([1.0; 1.0; 1.0], 0, t', lorenz)';

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
