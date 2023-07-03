function drawNodes(S ,nodes, timestep)

% collect nodes
xx = nodes.xx;
yy = nodes.yy;
mask_up = nodes.mask_up;
mask_down = nodes.mask_down;
mask_left = nodes.mask_left;
mask_right = nodes.mask_right;
Np = sum(sum(~isnan(mask_up)));
Nc = sum(sum(isnan(mask_left)));

% 
for k = 1:numel(S.h)
    if (k<=Np)
        % primary link
        sx = xx.*mask_up;
        sx(end,:) = [];
        sy = yy.*mask_up;
        sy(end,:) = [];
        sx = sx(:);
        sx(isnan(sx)) = [];
        sy = sy(:);
        sy(isnan(sy)) = [];
        
        ex = xx.*mask_down;
        ex(1,:) = [];
        ey = yy.*mask_down;
        ey(1,:) = [];
        ex = ex(:);
        ex(isnan(ex)) = [];
        ey = ey(:);
        ey(isnan(ey)) = [];
        set(S.h(k), 'XData', [sx(k),ex(k)], 'YData', [sy(k),ey(k)])
    else
        % cross-link
        sx = xx.*mask_right;
        sx(:,end) = [];
        sy = yy.*mask_right;
        sy(:,end) = [];
        sx = sx(:);
        sx(isnan(sx)) = [];
        sy = sy(:);
        sy(isnan(sy)) = [];
        
        ex = xx.*mask_left;
        ex(:,1) = [];
        ey = yy.*mask_left;
        ey(:,1) = [];
        ex = ex(:);
        ex(isnan(ex)) = [];
        ey = ey(:);
        ey(isnan(ey)) = [];
        set(S.h(k), 'XData', [sx(k-Np),ex(k-Np)], 'YData', [sy(k-Np),ey(k-Np)])
    end
end
set(S.mText,'String', timestep);
daspect([1,1,1])
drawnow;

end