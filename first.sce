// Brenton Horne's script 20161204
clear all

// Constants
// N for Chebyshev extrema grid
N         = 1000;
// NN for interpolation linspace grid
NN        = 10000;
// [a,b] is the problem interval
a         = 0;
b         = 200;
// - d2Y/dx2 + kxY = lambda Y is the problem; this is that k
k         = 1;
// How far we are computing the root-mean square error up to
M         = round(N/3);
format('v',16);

///////////////////////////////////////////////////////////////////////////////////////////////////
// Basic matrices                                                                                  //
///////////////////////////////////////////////////////////////////////////////////////////////////

// Vectors
/// n vector: a row vector of elements from 0 to N
n         = 0:N;
/// nsub vector: the n vector without 0 and N
nsub      = n(2:N);
/// Chebyshev extrema grid
x         = -cos(%pi*n'/N);
/// grid without endpoints
xsub      = x(2:N);
/// map Chebyshev extrema grid to integration interval: [a,b]
y         = (b-a)/2*x+(a+b)/2;
/// y without endpoints a and b
ysub      = (b-a)/2*xsub+(a+b)/2;
/// Linearly-spaced grid on [-1,1] for later interpolation
xx        = linspace(-1,1,NN+1)';
/// map interpolation grid to integration interval: [a,b]
yy        = (b-a)/2*xx+(a+b)/2;

///////////////////////////////////////////////////////////////////////////////////////////////////

// More square-shaped matrices
/// Chebyshev T(x) matrix for extrema grid
T         = cos(acos(x)*n);
/// T without x endpoints
Tsub      = T(2:N,:);
/// Chebyshev matrix for linearly-spaced interpolation grid
TT        = cos(acos(xx)*n);
/// Chebyshev U(x) matrix for extrema grid without endpoints (singularities exist there)
Usub      = diag(1./sqrt(1-xsub.^2))*sin(acos(xsub)*n);
/// First derivative of the Chebyshev T(x) function, on the extrema grid without endpoints
dTsub     = Usub*diag(n);
/// dT with endpoints
dT        = [-((-1).^n).*n.^2; dTsub; n.^2];
/// Second derivative of T(x) on extrema grid without endpoints
d2Tsub    = diag(1./(1-xsub.^2))*(diag(xsub)*Usub-Tsub*diag(n))*diag(n);
/// Second derivative of T(x) on extrema grid with endpoints
d2T       = [((-1).^n).*(n.^2).*(n.^2-1)/3; d2Tsub; (n.^2).*(n.^2-1)/3];
/// Second-order differentiation matrix for extrema grid
D2        = d2T/T;
/// Second-order differentiation matrix for extrema grid without endpoints
E2        = D2(2:N,2:N);

///////////////////////////////////////////////////////////////////////////////////////////////////
// Computation                                                                                     //
///////////////////////////////////////////////////////////////////////////////////////////////////

// H is the matrix version of the Sturm-Liouville operator.
// 4/(b-a)^2*E2 is the second-order differentiation matrix on the mapped interval
H         = - 4/((b-a)^2)*E2 + k*diag(ysub);
// Solution to specenvalue problem H Y = Lam Y
[Y, LAM]  = spec(H);
// Convert LAM (diagonal matrix) to vertical vector.
Lam       = diag(LAM);
// Order Lam in ascending order
[Lam, IX] = gsort(Lam, 'c', 'i');
// Order Y in the same order
Y         = Y(:,IX);
// Add boundary values
Y         = [zeros(1,N-1); Y; zeros(1,N-1)];

///////////////////////////////////////////////////////////////////////////////////////////////////
// Interpolate                                                                                     //
///////////////////////////////////////////////////////////////////////////////////////////////////

// Determine expansion coefficients for Chebyshev series
aa        = T\Y;

// Interpolate to linear grid
YY        = TT*aa;
// Amplitude
Amp       = YY.^2;

// Error analysis
// the exact solution calls for the specenvalues to negative the zeros of the Airy Ai function
//err       = airy(0,-Lam);
// root-mean square error for err vector up to element M
//rms       = sqrt(sum(err(1:M).^2)/(M));

///////////////////////////////////////////////////////////////////////////////////////////////////
// Plotting                                                                                        //
///////////////////////////////////////////////////////////////////////////////////////////////////

// Plot the log of err
//figure(1);
//plot(nsub,log10(err),"linewidth",2)
//title("Log10 of Airy Ai of - the specenvalues; lower (more negative) the better the solution")
// plot the first specenfunction
figure(1);
plot(yy,YY(:,1),'-r',"linewidth",2)
title("Plot of the first specenfunction on [0,200]")
// plot the 150th specenfunction
figure(2);
plot(yy,YY(:,150),'-g',"linewidth",2)
title("Plot of the 150th specenfunction on [0,200]")
// plot the 300th specenfunction
figure(3);
plot(yy,YY(:,300),'-m',"linewidth",2)
title("Plot of the 300th specenfunction on [0,200]")
// plot the 300th amplitude
figure(4)
plot(yy,Amp(:,300),'-k',"linewidth",2)
title("Plot of the 300th amplitude on [0,200]")
