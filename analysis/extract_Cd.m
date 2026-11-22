clear, clc, close all;
load(fullfile(fileparts(mfilename('fullpath')), '..', 'data', 'stage5_out.mat'));

t = t(1:end-1);
v = v(1:end-1);
m_grid = m_grid(1:end-1);

f = -D_hist ./ m_grid; % what the accel reads during coast. Up is positive
q = 0.5 * rho_hist .* v.^2;

k = 1:8:numel(f);

t = t(k);
f = f(k);
M_hist = M_hist(k);
q = q(k);
m_grid = m_grid(k);
Cd_hist = Cd_hist(k);

mask = (t > t_bo) & (M_hist >= 0.8) & (M_hist <= 1.2);
f_mask = f(mask); 
M_mask = M_hist(mask); 
q_mask = q(mask); 
m_mask = m_grid(mask); 
Cd_mask = Cd_hist(mask);

%% Least Squares
% checking unknowns
nodes = 0.80:0.05:1.20; % Row vector of Transonic Regime
dM = 0.05; % node spacing
n = numel(nodes); % number of unknowns
Ns = numel(f_mask); % number of samples in the masked force vector (transonic)
% creating empty system Ax = b
A = zeros(Ns, n); % creates a 1500 x 9 zero matrix
b = -f_mask(:); % f is negative because drag points downward. 
% filling A
for r = 1:Ns
    j = floor((M_mask(r) - nodes(1))/dM)+1; % interval where the Mach sample falls in. +1 because matlab counts from 1
    j = min(j, n-1); % caps j at 8 instead of 9, since 9+1 = 10, over the # of unknowns
    s = (M_mask(r) - nodes(j))/dM; % how far between the nodes do the sample sit
    wq = q_mask(r)*S/m_mask(r); % Shows how much deaccel one Cd unit produces at this sample's speed.
    A(r,j) = (1-s)*wq; % hat-function weight of left/right node, scaled from Cd to deceleration
    A(r,j+1) = s*wq; % same as comment above
end
% checking rows
for r = [1, 750, Ns] % checks first, middle, and last rows only
    fprintf('row %4d: nnz = %d | row sum = %.6f | qS/m = %.6f\n', ...
        r, nnz(A(r,:)), sum(A(r,:)), q_mask(r)*S/m_mask(r));
end
% solving two ways
c_qr = A\b; % Least squares 
c_ne = (A'*A)\(A'*b); % normal equations

% Diagnostics/Results
Cd_true = interp1(M_tab, Cd_tab, nodes, 'pchip');
err_pct = (c_qr' - Cd_true)./Cd_true*100;
residual = b - A*c_qr;
fprintf(['kappa(A) = %.4f | kappa(AtA) = ' ...
    '%.4f | kappa(A)^2 = %.4f\n'], ...
    cond(A), cond(A'*A), cond(A)^2);
fprintf('max |QR - NE| = %.3e\n', max(abs(c_qr - c_ne)));
fprintf('max |A''r| = %.3e\n', max(abs(A'*residual)));
fprintf('max |error| = %.4f %%\n', max(abs(err_pct)));

% Plotting
set(0, 'DefaultFigureVisible', 'on');
Mf = 0.80:0.001:1.20;
figure
plot(Mf, interp1(M_tab, Cd_tab, Mf, 'pchip'), 'k-', LineWidth = 1.5); hold on
plot(nodes, c_qr, 'ro-', LineWidth = 1)
xlabel('Mach'); ylabel('C_d'); grid on
legend('Truth (placeholder)', 'Recovered (least squares)', Location = 'northwest')
title('Cd(M) recovery, noise-free')