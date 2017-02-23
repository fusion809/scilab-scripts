// This is a script to solve the Airy SLE

N         = 1000;
M         = round(N/3);
NN        = 10000;
n         = 0:N;
nsub      = n(2:N);
nn        = 0:N;

// Initial conditions
a         = 0;
b         = 200;
k         = 1;

// Chebyshev grid
x         = cos(n'*%pi/N);
xx        = cos(nn'*%pi/N);
xsub      = x(2:N);

// map Chebyshev extrema grid to integration interval: [a,b]
y         = (b-a)/2*x+(a+b)/2;
// y without endpoints a and b
ysub      = (b-a)/2*xsub+(a+b)/2;
// Linearly-spaced grid on [-1,1] for later interpolation
xx        = linspace(-1,1,NN+1)';
// map interpolation grid to integration interval: [a,b]
yy        = (b-a)/2*xx+(a+b)/2;

// More square-shaped matrices
// Chebyshev T(x) matrix for extrema grid
T         = cos(acos(x)*n);
// T without x endpoints
Tsub      = T(2:N,:);

// T on the xx grid
TT        = cos(acos(xx)*n);
// Chebyshev U(x) matrix for extrema grid without endpoints (singularities exist there)
Usub      = diag(1./sqrt(1-xsub.^2))*sin(acos(xsub)*n);
// First derivative of the Chebyshev T(x) function, on the extrema grid without endpoints
dTsub     = Usub*diag(n);
// dT with endpoints
dT        = [-((-1).^n).*n.^2; dTsub; n.^2];
// Second derivative of T(x) on extrema grid without endpoints
d2Tsub    = diag(1./(1-xsub.^2))*(diag(xsub)*Usub-Tsub*diag(n))*diag(n);
// Second derivative of T(x) on extrema grid with endpoints
d2T       = [((-1).^n).*(n.^2).*(n.^2-1)/3; d2Tsub; (n.^2).*(n.^2-1)/3];
// Second-order differentiation matrix for extrema grid
D2        = d2T/T;
// Second-order differentiation matrix for extrema grid without endpoints
E2        = D2(2:N,2:N);

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// Computation                                                                                     //
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

//// H is the matrix version of the Sturm-Liouville operator.
//// 4/(b-a)^2*E2 is the second-order differentiation matrix on the mapped interval
H         = - 4/((b-a)^2)*E2 + k*diag(ysub);
//// Solution to eigenvalue problem H Y = Lam Y
[Y, LAM]  = spec(H);
//// Convert LAM (diagonal matrix) to vertical vector.
Lam       = diag(LAM);
//// Order Lam in ascending order
[Lam, IX] = gsort(Lam, 'r', 'i');
//// Order Y in the same order
Y         = Y(:,IX);
//// Add boundary values
Y         = [zeros(1,N-1); Y; zeros(1,N-1)];

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// Interpolate                                                                                     //
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

//// Determine expansion coefficients for Chebyshev series
aa        = T\Y;

//// Interpolate to linear grid
YY        = TT*aa;
//// Amplitude
Amp       = YY.^2;

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// Plotting                                                                                        //
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

figure(1)
plot(yy,Amp(:,100))
