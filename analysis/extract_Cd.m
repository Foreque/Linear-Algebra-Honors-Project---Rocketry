clear, clc, close all;
load(fullfile(fileparts(mfilename('fullpath')), '..', 'data', 'stage5_out.mat'));
set(0, 'DefaultFigureVisible', 'on'); % TURNS OFF/ON figures while coding

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
Mf = 0.80:0.001:1.20;
figure
plot(Mf, interp1(M_tab, Cd_tab, Mf, 'pchip'), 'w-', LineWidth = 1.5); hold on
plot(nodes, c_qr, 'ro-', LineWidth = 1)
xlabel('Mach'); ylabel('C_d'); grid on
legend('Truth (placeholder)', 'Recovered (least squares)', Location = 'northwest')
title('Cd(M) recovery, noise-free')

%% Noise, Weighted Least Squares, Covariance
% used ai to fetch information in ADXL357 datasheet, and also the equation
% ADXL357 datasheet Rev. A, Table 2, ±40 g
g0 = 9.80665; % force of gravity
nd_357 = 90e-6; % noise density at +_40gs
b_w_357 = 125; % bandwidth for Low Pass Filter Corner at 500 hz output data rate
lsb_357 = 78e-6; % Least Significant Bit (LSB); smallest measurable step the sensor can read
% total measurement noise standard deviation for an accelerometer
sigma_357 = g0 * sqrt((nd_357)^2*b_w_357+lsb_357^2/12);

% recreating the noisy measurements
rng(1); % fixes random seed, so I get the same "random" noise on every run
sig = sigma_357*ones(Ns,1); % 1500x1 column with one sigma per sample
b_n = b + sig .* randn(Ns, 1); % adding noise to b

% Weighted Least Squares
A_w  = A ./ sig; % divides every row by the standard deviation. Matlab stretches sig horizontally so that every row is divided
b_w  = b_n ./ sig; % divides every row to the noise vector b_n
c_w = A_w \ b_w; % solves for WLS using by solving linear least-squares via QR decomp
c_o = A \ b_n; % solves for unscaled OLS on the same data as comparison
% Computing the largets absolute difference between WLS coefficient and OLS
% coefficients, then print 
fprintf('max |WLS - OLS| = %.3e\n', max(abs(c_w - c_o)));

% Covariance and Error Bars
[~, R] = qr(A_w, "econ"); % Performing an QR decomp on weighted matrix A_w 
Rinv = R \ eye(n); % computing the inverse of the upper triangular matrix
Cov = Rinv*Rinv'; % solves the covariance matrix
se = sqrt(diag(Cov)); % standard error; takes the square root of the diagonal of the covariance matrix
                      % diagonal entries represent the variance of each
                      % individual parameter
% Computes the relative difference between the QR-derived covariance matrix
% and the textbook inverse using matrix norms
fprintf('rel. diff Cov vs inv(A_w''A_w) = %.3e\n', ...
    norm(Cov - inv(A_w'*A_w))/norm(Cov)); 

% Diagnostics 
% c_qr is the noise free solve. Subtracting it from the WLS shows how much
% error the noise caused
err_noise = c_w - c_qr; 
% Divides each paramaters noise error by its estimated standard error
z = err_noise ./ se;
% Computes Weighted Residuals, differnece between measured data and the models
% predictions
r_w = b_w - A_w*c_w;
% Average squared weighted residual. 
chi2_nu = (r_w'*r_w)/(Ns - n);

fprintf('max se/Cd = %.4f %%\n', max(se./Cd_true(:))*100);
fprintf('max |z| = %.2f\n', max(abs(z)));
fprintf('chi2/nu = %.1f\n', chi2_nu);
fprintf('max |total error| = %.4f %%\n', max(abs((c_w' - Cd_true)./Cd_true))*100);

% Plot
figure
errorbar(nodes(:), 100*err_noise./Cd_true(:), 100*1.96*se./Cd_true(:), 'wo'); 
yline(0, 'k--');
xlabel('Mach'); 
ylabel('Noise-only error (% of C_d)'); 
title('Noise error vs. predicted 95% interval (ADXL357)')
hold on
grid on