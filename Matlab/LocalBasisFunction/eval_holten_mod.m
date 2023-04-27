% Main property evaluation function
function results=eval_holten_mod(PT) % T in K, P in MPa

if(iscell(PT))
    P=PT{1};
    np=length(P);
    T=PT{2};
    nt=length(T);
    [pp,tt]=ndgrid(P,T);
    P=pp(:);
    T=tt(:);
    flgcell=1;
    npm=length(P);
else
    flgcell=0;
    P=PT(:,1);
    T=PT(:,2);
    npm=length(P);
end
P=P*1e6;

G=zeros(npm,1);
rho=G;
Kap=G;
CP=G;
U=G;
K=G;

% Parameters
Tc = 228.2; Pc = 0;
rho0 = 1081.6482;
R = 461.523087;
omega0 = 0.52122690;
L0 = 0.76317954;
k0 = 0.072158686;
k1 = -0.31569232;
k2 = 5.2992608;

% Background coefficients
c = [-8.1570681381655 1.2875032e+000 7.0901673598012 ...
	-3.2779161e-002 7.3703949e-001 -2.1628622e-001 -5.1782479e+000 ...
	4.2293517e-004 2.3592109e-002 4.3773754e+000 -2.9967770e-003 ...
	-9.6558018e-001 3.7595286e+000 1.2632441e+000 2.8542697e-001 ...
	-8.5994947e-001 -3.2916153e-001 9.0019616e-002 8.1149726e-002 ...
	-3.2788213e+000];
a = [0 0 1 -0.2555 1.5762 1.64 3.6385 -0.3828 ...
	1.6219 4.3287 3.4763 5.1556 -0.3593 5.0361 2.9786 6.2373 ...
	4.046 5.3558 9.0157 1.2194];
b = [0 1 0 2.1051 1.1422 0.951 0 3.6402 ...
	2.076 -0.0016 2.2769 0.0008 0.3706 -0.3975 2.973 -0.318 ...
	2.9805 2.9265 0.4456 0.1298];
d = [0 0 0 -0.0016 0.6894 0.013 0.0002 0.0435 ...
	0.05 0.0004 0.0528 0.0147 0.8584 0.9924 1.0041 1.0961 ...
	1.0228 1.0303 1.618 0.5213];

% Background derivatives

function res = B(tau, pi)
    res = sum(c .* tau.^a .* pi.^b .* exp(-d.*pi));
end

function res = Bp(tau, pi)
    res = sum(c .* tau.^a .* pi.^(b-1) .* (b - d .* pi) .* exp(-d .* pi));
end

function res = Bt(tau, pi)
    res = sum(c .* a .* tau.^(a-1) .* pi.^b .* exp(-d .* pi));
end

function res = Bpp(tau, pi)
    res = sum(c .* tau.^a .* pi.^(b-2) .* ((d .* pi - b).^2 - b) .* exp(-d .* pi));
end

function res = Btp(tau, pi)
    res = sum(c .* a .* tau.^(a-1) .* pi.^(b-1) .* (b - d .* pi) .* exp(-d .* pi));
end

function res = Btt(tau, pi)
    res = sum(c .* a .* (a-1) .* tau.^(a-2) .* pi.^b .* exp(-d .* pi));
end

P0 = -300e6;
G=zeros(npm,1);
S=G;
CP=G;
CV=G;
rho=G;
Kap=G;
Alp=G;
U=G;
E=G;

for i=1:npm
% Dimensionless temperature and pressure
t =  (T(i) - Tc) / Tc;
p =  (P(i) - Pc) / (rho0 * R * Tc);
tau = T(i) / Tc;
pi = (P(i) - P0) / (rho0 * R * Tc);

% Field L and its derivatives
K1 = sqrt((1+k0*k2+k1*(p-k2*t))^2 - 4*k0*k1*k2*(p-k2*t));
K3 = K1^3;
K2 = sqrt(1 + k2^2);
L = L0 * K2 * (1 - K1 + k0*k2 + k1*(p + k2*t)) / (2*k1*k2);
Lt = L0 * 0.5 * K2 * (1 + (1 - k0*k2 + k1*(p - k2*t))/K1);
Lp = L0 * K2 * (K1 + k0*k2 - k1*p + k1*k2*t - 1) / (2*k2*K1);
Ltt = -2*L0*K2*k0*k1*k2^2 / K3;
Ltp =  2*L0*K2*k0*k1*k2 / K3;
Lpp = -2*L0*K2*k0*k1 / K3;

% Interaction parameter omega
omega = 2 + omega0 * p;

% Calculate equilibrium fraction xe
x = findxe(L, omega);

% Order parameter f and susceptibility chi
f = 2 * x - 1;
f2 = f^2;
chi = 1 / (2 / (1 - f2) - omega);

% Dimensionless properties
g0 = x*L + x*log(x) + (1-x)*log(1-x) + omega*x*(1-x); % g0 = (g - gA)/tau
g = B(tau,pi) + tau*g0;
s = -0.5*(f+1)*Lt*tau - g0 - Bt(tau,pi);
v = 0.5 * tau * (omega0/2*(1-f2) + Lp*(f+1)) + Bp(tau,pi);
kap = (1/v) * (tau/2 * (chi * (Lp - omega0 * f)^2 - (f+1)*Lpp) - Bpp(tau,pi));
alp = (1/v) * (Ltp/2 * tau*(f+1) + (omega0/2*(1-f2) + Lp*(f+1))/2 ...
    - tau*Lt/2 * chi*(Lp - omega0*f) + Btp(tau,pi));
cp = tau * ( -Lt * (f+1) + tau*(Lt^2 * chi - Ltt*(f+1)) / 2 - Btt(tau,pi));

% Properties in SI units
G(i) = R * Tc * g;                     % Specific Gibbs energy
S(i) = R * s;							% Specific entropy
rho(i) = rho0 / v;						% Density
Kap(i) = kap / (rho0 * R * Tc);		% Isothermal compressibility
Alp(i) = alp / Tc;						% Expansion coefficient
CP(i) = R * cp;						% Isobaric heat capacity
CV(i) = CP(i) - T(i)*Alp(i)^2 / (rho(i) * Kap(i));	% Isochoric heat capacity
U(i) = 1/sqrt(rho(i)*Kap(i) - T(i)*Alp(i)^2/CP(i));% Speed of sound
K(i)=1e-9/Kap(i);
E(i)=G(i)-P(i)/rho(i)+T(i)*S(i);
end

if(flgcell)
    G=reshape(G,np,nt);
    rho=reshape(rho,np,nt);
    K=reshape(K,np,nt);
    CP=reshape(CP,np,nt);
    CV=reshape(CV,np,nt);
    S=reshape(S,np,nt);
    U=reshape(U,np,nt);
    E=reshape(E,np,nt);
    Alp=reshape(Alp,np,nt);
end

% 
% fprintf(['Temperature   \t%f\tK\n' ...
%          'Pressure      \t%f\tMPa\n' ...
%          'Fraction xe   \t%.8f\t\n' ...
%          'Field L       \t%.8f\t\n' ...
%          'Gibbs energy  \t%f\tJ/kg\n' ...
%          'Entropy       \t%f\tJ/(kg K)\n' ...
%          'Density       \t%.5f\tkg/m3\n' ...
%          'Isoth. comp.  \t%e\t1/MPa\n' ...
%          'Expansivity   \t%e\t1/K\n' ...
%          'CP            \t%.4f\tJ/(kg K)\n' ...
%          'CV            \t%.4f\tJ/(kg K)\n' ...
%          'Speed of sound\t%.4f\tm/s\n\n'], ...
%          [T, 1e-6*P, x, L, G, S, rho, 1e6 * Kap, Alp, CP, CV, U]);

results.rho=rho;
results.G=G;
results.K=K;
results.Cp=CP;
results.Cv=CV;
results.S=S;
results.E=E;
results.vel=U;
results.alpha=Alp;
end

% Function for calculating the equilibrium fraction xe
function xe = findxe(L, W)
    flip =  L < 0;
    if flip
        L = -L;
    end
    
    % Find starting values for x
    if W < 1.1111111*(2.944439 - L) % xe = 0.05 isoline, W = (10/9) * (ln(19) - L)
        x0 = 0.049;
        x1 = 0.5;
    elseif W < 1.0204081*(4.595119 - L) % xe = 0.01 isoline, W = (50/49) * (ln(99) - L)
        x0 = 0.0099;
        x1 = 0.051;
    else
        x0 = 0.99 * exp(-1.0204081 * L - W);
        x1 = 1.01 * 1.087 * exp(-L - W);
        if x1 > 0.0101
            x1 = 0.0101;
        end
    end

    xefun = @(x) L + log(x/(1-x)) + W*(1-2*x);
    
    xe = fzero(xefun, [x0 x1]);
    if flip
        xe = 1 - xe;
    end
    
end
    