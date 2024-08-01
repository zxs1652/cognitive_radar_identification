function I_xy = mutualInfo(x,y,pts1,pts2)
%MUTUALINFO This function is used to numerically estimate the mutual
%information between two variable x, y
%   Input  -- x, y: data points of two variables
%             pts1, pts2: resampling points of the distribution
%   Output -- I_xy: mutual information of x, y

% estimate marginal distribution (KDE)
x_min = min(x); x_max = max(x);
% if abs(x_min - x_max) < 4.5e-13
%     pts1 = x_min;
% else
%     pts1 = linspace(x_min, x_max, 200);
% end
pts1 = linspace(x_min, x_max, 200);

[density_px, xi_px] = ksdensity(x, pts1);
density_px = density_px / trapz(xi_px, density_px) * 1;
density_px(density_px < 1e-20) = 1e-20;
[density_py, xi_py] = ksdensity(y, pts2);
density_py = density_py / trapz(xi_py, density_py) * 1;
density_py(density_py < 1e-20) = 1e-20;

% estimate joint distribution (KDE)
[grid1,grid2] = meshgrid(pts1, pts2);
grid_pdf = [grid1(:), grid2(:)];    
[density, ~] = ksdensity([x; y]', grid_pdf);
density_joint = density';

% compute numerical double integral for KL divergence
density_joint_2d = reshape(density_joint, length(pts2), length(pts1));
density_joint_2d = density_joint_2d / trapz(xi_px,trapz(xi_py,density_joint_2d,1)) * 1;
density_joint_2d(density_joint_2d < 1e-40) = 1e-40;
px_2d = repmat(density_px, length(pts2), 1);
py_2d = repmat(density_py', 1, length(pts1));
f_mi = density_joint_2d.*(log(density_joint_2d./(px_2d.*py_2d)));
f_mi(isnan(f_mi)) = 0;
I_xy= trapz(xi_py,trapz(xi_px,f_mi,2));


end

