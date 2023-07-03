function nodes = updateNode(nodes, mass, l0, stiffness, stiffness2, damping, w, dt)
r0 = w/2;                                                                   % the cellular radius

% collect nodes
xx = nodes.xx;
yy = nodes.yy;
vel_xx = nodes.vel_xx;
vel_yy = nodes.vel_yy;
isFixed = nodes.isFixed;
mask_up = nodes.mask_up;
mask_down = nodes.mask_down;
mask_left = nodes.mask_left;
mask_right = nodes.mask_right;

% calculate force on each node
% up
dx = zeros(size(xx));
dx(1:end-1,:) = -diff(xx,1,1);
dy = zeros(size(yy));
dy(1:end-1,:) = -diff(yy,1,1);
norm_dxdy = sqrt(dx.^2+dy.^2);
dv_x = zeros(size(vel_xx));
dv_x(1:end-1,:) = -diff(vel_xx,1,1);
dv_y = zeros(size(vel_yy));
dv_y(1:end-1,:) = -diff(vel_yy,1,1);
fx_up = (-stiffness2 * (norm_dxdy - l0) .* (dx ./ norm_dxdy) - (damping*dv_x)).*mask_up;
fx_up(isnan(fx_up)) = 0;
fy_up = (-stiffness2 * (norm_dxdy - l0) .* (dy ./ norm_dxdy) - (damping*dv_y)).*mask_up;
fy_up(isnan(fy_up)) = 0;

% down
dx = zeros(size(xx));
dx(2:end,:) = diff(xx,1,1);
dy = zeros(size(yy));
dy(2:end,:) = diff(yy,1,1);
norm_dxdy = sqrt(dx.^2+dy.^2);
dv_x = zeros(size(vel_xx));
dv_x(2:end,:) = diff(vel_xx,1,1);
dv_y = zeros(size(vel_yy));
dv_y(2:end,:) = diff(vel_yy,1,1);
fx_down = (-stiffness2 * (norm_dxdy - l0) .* (dx ./ norm_dxdy) -damping*dv_x).*mask_down;
fx_down(isnan(fx_down)) = 0;
fy_down = (-stiffness2 * (norm_dxdy - l0) .* (dy ./ norm_dxdy) -damping*dv_y).*mask_down;
fy_down(isnan(fy_down)) = 0;

% left
dx = zeros(size(xx));
dx(:,2:end) = diff(xx,1,2);
dy = zeros(size(yy));
dy(:,2:end) = diff(yy,1,2);
norm_dxdy = sqrt(dx.^2+dy.^2);
dv_x = zeros(size(vel_xx));
dv_x(:,2:end) = diff(vel_xx,1,2);
dv_y = zeros(size(vel_yy));
dv_y(:,2:end) = diff(vel_yy,1,2);
fx_left = (-stiffness * (norm_dxdy - l0) .* (dx ./ norm_dxdy) + (-damping*dv_x)).*mask_left;
fx_left(isnan(fx_left)) = 0;
fy_left = (-stiffness * (norm_dxdy - l0) .* (dy ./ norm_dxdy)+ (-damping*dv_y)).*mask_left; 
fy_left(isnan(fy_left)) = 0;
         
% right
dx = zeros(size(xx));
dx(:,1:end-1) = -diff(xx,1,2);
dy = zeros(size(yy));
dy(:,1:end-1) = -diff(yy,1,2);
norm_dxdy = sqrt(dx.^2+dy.^2);
dv_x = zeros(size(vel_xx));
dv_x(:,1:end-1) = -diff(vel_xx,1,2);
dv_y = zeros(size(vel_yy));
dv_y(:,1:end-1) = -diff(vel_yy,1,2);
fx_right = (-stiffness * (norm_dxdy - l0) .* (dx ./ norm_dxdy) + (-damping*dv_x)).*mask_right;
fx_right(isnan(fx_right)) = 0;
fy_right = (-stiffness * (norm_dxdy - l0) .* (dy ./ norm_dxdy)+ (-damping*dv_y)).*mask_right;
fy_right(isnan(fy_right)) = 0;

% total force
force_xx =  fx_up + fx_down + fx_left + fx_right;
force_yy =  fy_up + fy_down + fy_left + fy_right;

% update
acc_xx = (force_xx/mass).*(~isFixed);
vel_xx = (vel_xx+acc_xx*dt).*(~isFixed);
xx = xx+vel_xx*dt;
acc_yy = (force_yy/mass).*(~isFixed);
vel_yy = (vel_yy+acc_yy*dt).*(~isFixed);
yy = yy+vel_yy*dt;

% virtual term correction
force_xx(end,:) = force_xx(2,:);
force_xx(1,:) = force_xx(end-1,:);
force_yy(end,:) = force_yy(2,:);
force_yy(1,:) = force_yy(end-1,:);
acc_xx(end,:) = acc_xx(2,:);
acc_xx(1,:) = acc_xx(end-1,:);
acc_yy(end,:) = acc_yy(2,:);
acc_yy(1,:) = acc_yy(end-1,:);
vel_xx(end,:) = vel_xx(2,:);
vel_xx(1,:) = vel_xx(end-1,:);
vel_yy(end,:) = vel_yy(2,:);
vel_yy(1,:) = vel_yy(end-1,:);
xx(end,:) = xx(2,:);
xx(1,:) = xx(end-1,:);
yy(end,:) = yy(2,:)+2*pi*r0;
yy(1,:) = yy(end-1,:)-2*pi*r0;

% renew nodes
nodes.xx = xx;
nodes.yy = yy;
nodes.vel_xx = vel_xx;
nodes.vel_yy = vel_yy;
nodes.acc_xx = acc_xx;
nodes.acc_yy = acc_yy;
nodes.force_xx = force_xx;
nodes.force_yy = force_yy;

