clear, clc, close all;
load(fullfile(fileparts(mfilename('fullpath')), '..', 'data', 'stage5_out.mat'));

t = t(1:end-1);
v = v(1:end-1);
m_grid = m_grid(1:end-1);

f = -D_hist ./ m_grid;
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