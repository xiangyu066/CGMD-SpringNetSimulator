function nodes = buildNodes(w, L, rho_lattice, rho_loop)

% conver width and length to cylinder and semi-sphere
r0 = w/2;                                                                   % the cellular radius
h = L - 2*r0;                                                               % height of cylinder part

% cylinder part
nLoops_cylinder = round(rho_loop*h);            
nLattices_cylinder = round(rho_lattice*2*pi*r0);

% create meshgrid
rows = nLattices_cylinder;
cols = nLoops_cylinder;
[cc,rr] = meshgrid(1:cols,1:rows);
xx0 = (h/(cols-1))*(cc-1);
yy0 =  (2*pi*r0/(rows-1))*(rr-1);

% initial condition
xx0(end,:) = [];                                                            % remove overlap terms between beginning and final term because periodic boundary
xx0 = [xx0(end,:);xx0;xx0(1,:)];                                            % add virtual terms
yy0(end,:) = [];                                                            % remove overlap terms between beginning and final term because periodic boundary
yy0 = [yy0(end,:)-2*pi*r0;yy0;yy0(1,:)+2*pi*r0];                            % add virtual terms
xx = xx0;
yy = yy0;
vel_xx = zeros(size(xx));
vel_yy = zeros(size(yy));
acc_xx = zeros(size(xx));
acc_yy = zeros(size(yy));
force_xx = zeros(size(xx));
force_yy = zeros(size(yy));
isFixed = ones(size(xx));
isFixed(:,2:end-1) = 0;

% check even or odd number
rows = rows-1;                                                              % remove overlap terms between beginning and final term because periodic boundary
if mod(rows,2)>0
    rows_ = rows+1;
else
    rows_ = rows;
end

if mod(cols,2)>0
    cols_ = cols+1;
else
    cols_ = cols;
end

% label the left-hand side grid
kernal = [NaN,1;1,NaN];
mask_left = repmat(kernal,[rows_,cols_]/2);
mask_left(:,1) = 0;                                                         % correct after consider semi sphere

if rows_~=rows
   mask_left(end,:) = []; 
end

if cols_~=cols
   mask_left(:,end) = []; 
end

% label the right-hand side grid
kernal = [1,NaN;NaN,1];
mask_right = repmat(kernal,[rows_,cols_]/2);
mask_right(:,end) = 0;                                                      % correct after consider semi sphere

if rows_~=rows
   mask_right(end,:) = []; 
end

if cols_~=cols
   mask_right(:,end) = []; 
end

% label primary link
mask_up = nan(rows,cols);
for col = 1:cols
    row =  randi([1,rows],1);
    idxs = (1:rows);                                                        % create row listing
    idxs = circshift(idxs,-(row-1));                                        % shift beginning row
    while (~isempty(idxs))
        strands_num = randi([5,9],1);                                       % determine lengths of glycan strands
        if (strands_num<length(idxs))
            mask_up(idxs(1:strands_num),col) = 1;
            idxs(1:strands_num+1) = [];
        else
            if (length(idxs)>1)
                mask_up(idxs(1:end-1),col) = 1;
                idxs(1:end) = [];
            else
                isgroup = randi([0,1],1);
                if (isgroup==1)
                    mask_up(idxs-1,col) = 1;
                else
                    mask_up(idxs,col) = 1;
                end
                idxs = [];
            end
        end
    end
end
mask_down = circshift(mask_up,[1,0]);

% add virtual terms
mask_up = [mask_up(end,:);mask_up;mask_up(1,:)];
mask_up(end,:) = NaN;
mask_down = [mask_down(end,:);mask_down;mask_down(1,:)];
mask_down(1,:) = NaN;
mask_left = [mask_left(end,:);mask_left;mask_left(1,:)];
mask_right = [mask_right(end,:);mask_right;mask_right(1,:)];

% construct nodes
nodes.xx0 = xx0;
nodes.yy0 = yy0;
nodes.xx = xx;
nodes.yy = yy;
nodes.vel_xx = vel_xx;
nodes.vel_yy = vel_yy;
nodes.acc_xx = acc_xx;
nodes.acc_yy = acc_yy;
nodes.force_xx = force_xx;
nodes.force_yy = force_yy;
nodes.isFixed = isFixed;
nodes.mask_up = mask_up;
nodes.mask_down = mask_down;
nodes.mask_left = mask_left;
nodes.mask_right = mask_right;

end
