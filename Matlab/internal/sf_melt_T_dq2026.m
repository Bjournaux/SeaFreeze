function Tm = sf_melt_T_dq2026(P)
% SF_MELT_T_DQ2026  Melting temperature (K) of the stable solid at P (MPa).
%
% Port of psiEOS.m melt_T (lbf-thermo, model 'dq2026', JMB 2026): IAPWS
% R14-08 for ice Ih/III/V/VI to 2.17 GPa; ice VII from Datchi et al. (2000)
% and superionic VII'' from Queyroux et al. (2020) to 45 GPa; above, linear
% in ln P onto French & Hamel's superionic -> fluid line.  273.16 K below the
% triple-point pressure.  Used as the water_Brown2026 validity mask in sf_phase_map
% (the surface is not water more than 40 K below this curve).
sz = size(P); P = P(:); Tm = NaN(size(P));
PT_si = [45 850 * ((45 - 14.6) / 3.44 + 1)^(1 / 4.33); 330.6 6000; 627.6 7000; 1034.2 8000; ...
         1592.9 10000; 3734.6 10000; 4664.8 12000; 9311.3 12000];
Tm(P < 611.657e-6) = 273.16;
ih = P >= 611.657e-6 & P < 208.566;
a = 251.165 * ones(nnz(ih), 1); b = 273.16 * ones(nnz(ih), 1); Pi = P(ih);
for k = 1:60, m = 0.5 * (a + b); up = p_melt_Ih(m) > Pi; a(up) = m(up); b(~up) = m(~up); end
Tm(ih) = 0.5 * (a + b);
hp = P >= 208.566 & P < 2170 & isfinite(P);
a = 251.165 * ones(nnz(hp), 1); b = 2e4 * ones(nnz(hp), 1); Ph = P(hp);
for k = 1:80, m = 0.5 * (a + b); up = p_melt_hp(m) < Ph; a(up) = m(up); b(~up) = m(~up); end
Tm(hp) = 0.5 * (a + b);
vii = P >= 2170 & P <= 45000 & isfinite(P);
Pg = P(vii) / 1e3;
Td = 354.8 * max((Pg - 2.17) / 1.253 + 1, 1e-12).^(1/3);
xq = (Pg - 14.6) / 3.44 + 1; Tq = zeros(size(Pg)); Tq(xq > 0) = 850 * xq(xq > 0).^(1/4.33);
Tm(vii) = max(Td, Tq);
si = P > 45000 & isfinite(P);
Tm(si) = interp1(log(PT_si(:, 1)), PT_si(:, 2), log(min(P(si) / 1e3, PT_si(end, 1))), 'linear');
Tm = reshape(Tm, sz);
end

function p = p_melt_Ih(T)
th = T / 273.16;
p = 611.657e-6 * (1 + 0.119539337e7 * (1 - th.^3) + 0.808183159e5 * (1 - th.^25.75) + 0.333826860e4 * (1 - th.^103.75));
end

function p = p_melt_hp(T)
p = NaN(size(T));
m = T >= 251.165 & T < 256.164; th = T(m) / 251.165; p(m) = 208.566 * (1 - 0.299948 * (1 - th.^60));
m = T >= 256.164 & T < 273.31;  th = T(m) / 256.164; p(m) = 350.1 * (1 - 1.18721 * (1 - th.^8));
m = T >= 273.31 & T < 355.0;    th = T(m) / 273.31;  p(m) = 632.4 * (1 - 1.07476 * (1 - th.^4.6));
m = T >= 355.0 & T <= 715.0;    th = T(m) / 355.0;
p(m) = 2216.0 * exp(0.173683e1 * (1 - 1 ./ th) - 0.544606e-1 * (1 - th.^5) + 0.806106e-7 * (1 - th.^22));
t7 = 715 / 355; P715 = 2216.0 * exp(0.173683e1 * (1 - 1 / t7) - 0.544606e-1 * (1 - t7^5) + 0.806106e-7 * (1 - t7^22));
n = log(45000 / P715) / log(1600 / 715);
m = T > 715.0; p(m) = P715 * (T(m) / 715).^n;
end
