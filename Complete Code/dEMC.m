function [MieC, iImage] = dEMC(MieScatter, Light, Dynamic, G, Camera)
% Monte Carlo simulation for polarized light transport in scattering media
% Inputs:
%   MieScatter - Structure containing Mie scattering parameters
%   Light      - Structure containing light source parameters
%   Dynamic    - Structure containing dynamic properties
%   G          - Structure containing geometry information
%   Camera     - Structure containing camera parameters
% Outputs:
%   MieC       - Structure containing Mie scattering coefficients
%   iImage     - Structure containing simulation results (16 variables)

format long;

% ==================== 1. Extract parameters from input structures ====================

% Extract parameters from Light structure
PhotonNum = Light.PhotonNum;
c0 = Light.c0;
lambda = Light.lambda;

% Extract parameters from Camera structure
snap = Camera.snap;
timeScale = snap * 1e6; % Convert to microseconds
Size_per_Pixel = Camera.xSize;
zSize = Camera.zSize;
pol_angles = Camera.pol_angles;
num_pol_angles = length(pol_angles);
NA = Camera.NA;

% Extract parameters from MieScatter structure
n_water = MieScatter.n_water;
NN = MieScatter.NN;
MieTypeNum = MieScatter.MieTypeNum;
n_sctCOMPLEX = MieScatter.n_sctCOMPLEX;
ua = MieScatter.ua;
dia = MieScatter.dia;
Volfrac = MieScatter.Volfrac;

% Precompute Mie scattering parameters
for ll = 1:MieTypeNum
    [xSize(ll), qext(ll), qsca(ll), albedo(ll), gFactor(ll), psca0(ll), s1tab(ll,:), s2tab(ll,:)] = ...
        PreparMieScattering(lambda, n_water, n_sctCOMPLEX(ll), dia(ll), NN);
end

% Calculate optical properties
Lparticle = (4*pi./(3*Volfrac)).^(1/3).*(dia/2);
lsa = 2*dia./(3.*Volfrac.*qext);
Klsa = 2*pi*n_water/lambda;
lt = lsa./(1 - gFactor);

% Extract parameters from Geometry structure
X = G.X; % Number of voxels (x-direction)
Y = G.Y; % Number of voxels (y-direction)
L = G.L; % Number of voxels (z-direction)
MediaTypeMapping = G.MediaTypeMapping;

% Determine media complexity level
if isfield(G, 'complex')
    complexity = G.complex;
else
    % Automatically determine complexity based on scattering mean free path and voxel size
    lt_min = min(lt);
    voxel_size = Camera.Size_per_Pixel;
    ratio = lt_min / voxel_size;
    
    if ratio > 100
        complexity = 0; % Homogeneous media
    elseif ratio > 10
        complexity = 1; % Simple boundary detection
    elseif ratio > 3
        complexity = 2; % Medium boundary detection
    else
        complexity = 3; % Complex boundary detection
    end
end

fprintf('Media complexity level: %d\n', complexity);

% Extract parameters from Dynamic structure
VelocityMappingx = Dynamic.FVxMapping;
VelocityMappingy = Dynamic.FVyMapping;
VelocityMappingz = Dynamic.FVzMapping;
D = Dynamic.BMMapping;

% Generate random directions for Brownian motion
rand_dir_x_1 = 2*rand(X, Y, L)-1;
rand_dir_y_1 = 2*rand(X, Y, L)-1;
rand_dir_z_1 = 2*rand(X, Y, L)-1;

% Normalize random directions
norm_factor = sqrt(rand_dir_x_1.^2 + rand_dir_y_1.^2 + rand_dir_z_1.^2);
rand_dir_x_1 = rand_dir_x_1 ./ norm_factor;
rand_dir_y_1 = rand_dir_y_1 ./ norm_factor;
rand_dir_z_1 = rand_dir_z_1 ./ norm_factor;

% Calculate Brownian motion velocity components (based on Einstein relation)
v_brownian = sqrt(2*D/timeScale); % Brownian motion characteristic velocity
VelocityMappingx = VelocityMappingx + v_brownian .* rand_dir_x_1;
VelocityMappingy = VelocityMappingy + v_brownian .* rand_dir_y_1;
VelocityMappingz = VelocityMappingz + v_brownian .* rand_dir_z_1;

% Calculate maximum acceptance angle for detection
max_angle = asin(NA / n_water); % Maximum acceptance angle (radians)
cos_max_angle = cos(max_angle); % For fast comparison

% Store optical properties in MieC structure
MieC.n_sctCOMPLEX = n_sctCOMPLEX;
MieC.dia = dia;
MieC.Volfrac = Volfrac;
MieC.lsa = lsa;
MieC.lt = lt;
MieC.gFactor = gFactor;
MieC.complexity = complexity; % Store complexity level

% ==================== 2. Initialize variables ====================

% Photon state variables
MAXSCT = 1000;
dd = -log(rand(PhotonNum, MAXSCT));

% Calculate image dimensions
scaled_X = X;
scaled_Y = Y;

% Initialize the 16 output variables (8 for reflection, 8 for transmission)

% Reflection variables (3D, with layers)
IMRS_layer_pol = zeros(scaled_X, scaled_Y, L, num_pol_angles);  % Single scattering reflection polarized images
IMRM_layer_pol = zeros(scaled_X, scaled_Y, L, num_pol_angles);  % Multiple scattering reflection polarized images
IMRS_layer = zeros(scaled_X, scaled_Y, L);                      % Single scattering reflection non-polarized images
IMRM_layer = zeros(scaled_X, scaled_Y, L);                      % Multiple scattering reflection non-polarized images
G1_perm_RS_layer = zeros(scaled_X, scaled_Y, L);               % Single scattering reflection g1 values
G1_perm_RM_layer = zeros(scaled_X, scaled_Y, L);               % Multiple scattering reflection g1 values
Photonpath_perm_RS_layer = zeros(scaled_X, scaled_Y, L);       % Single scattering reflection photon path counts
Photonpath_perm_RM_layer = zeros(scaled_X, scaled_Y, L);       % Multiple scattering reflection photon path counts

% Transmission variables (3D, with layers)
IMTS_layer_pol = zeros(scaled_X, scaled_Y, L, num_pol_angles);  % Single scattering transmission polarized images
IMTM_layer_pol = zeros(scaled_X, scaled_Y, L, num_pol_angles);  % Multiple scattering transmission polarized images
IMTS_layer = zeros(scaled_X, scaled_Y, L);                      % Single scattering transmission non-polarized images
IMTM_layer = zeros(scaled_X, scaled_Y, L);                      % Multiple scattering transmission non-polarized images
G1_perm_TS_layer = zeros(scaled_X, scaled_Y, L);               % Single scattering transmission g1 values
G1_perm_TM_layer = zeros(scaled_X, scaled_Y, L);               % Multiple scattering transmission g1 values
Photonpath_perm_TS_layer = zeros(scaled_X, scaled_Y, L);       % Single scattering transmission photon path counts
Photonpath_perm_TM_layer = zeros(scaled_X, scaled_Y, L);       % Multiple scattering transmission photon path counts

% ==================== 3. Main photon transport loop ====================
disp('Begin Photon Transport Cycle');
tic;

for kk = 1:PhotonNum
    % Display progress
    if mod(kk, 10000) == 0
        fprintf('Processing photon %d/%d\n', kk, PhotonNum);
    end
    
    % Initialize photon
    Ei1 = complex(1, 0);
    Ei2 = 0;
    l = 1; m = 0; n = 0;
    AngleSita = 0;
    AnglePhi = 0;
    u = sin(AngleSita)*cos(AnglePhi);
    v = sin(AngleSita)*sin(AnglePhi);
    w = cos(AngleSita);
    p = n*v - m*w;
    q = l*w - n*u;
    r = m*u - l*v;
    
    % Initialize polarization matrix
    P = zeros(3,3);
    P(1,1) = l; P(1,2) = m; P(1,3) = n;
    P(2,1) = p; P(2,2) = q; P(2,3) = r;
    P(3,1) = u; P(3,2) = v; P(3,3) = w;
    
    nsct = 1;
    
    % Initialize photon state variables
    t = 0;
    rejmu = 0;
    rejphi = 0;
    
    % Initial position (random in x-y plane, z=0)
    x = (rand(1)) * X * Size_per_Pixel - 0.5 * X * Size_per_Pixel;
    y = (rand(1)) * Y * Size_per_Pixel - 0.5 * Y * Size_per_Pixel;
    z = 0;
    s = 0;
    g1 = 1; % Initial g1 value
    dphase = 0;
    
    weight = 1;
    E1 = Ei1;
    E2 = Ei2;
    
    % Photon transport loop
    while true
        % Check if maximum scattering events reached
        if nsct >= MAXSCT
            break;
        end
        
        % Store current position for boundary detection
        x_prev = x;
        y_prev = y;
        z_prev = z;
        t_prev = t;
        
        % Get current tissue type
        [indx, indy, indz] = tissueTypeMapping(x/Size_per_Pixel, y/Size_per_Pixel, z/zSize, X, Y, L);
        current_tissueType = MediaTypeMapping(indx, indy, indz);
        
        % Boundary detection and movement
        dd_current = dd(kk, nsct);
        
        switch complexity
            case 0
                % Homogeneous media, no boundary detection
                [x_new, y_new, z_new, crossed_boundary, boundary_info] = movePhotonNoBoundaryCheck(...
                    x, y, z, P, dd_current, lsa(current_tissueType), current_tissueType);
                
            case 1
                % Simple boundary detection - bisection method
                [x_new, y_new, z_new, crossed_boundary, boundary_info] = movePhotonSimpleBoundary(...
                    x, y, z, P, dd_current, lsa, MediaTypeMapping, Size_per_Pixel, X, Y, L, current_tissueType);
                
            case 2
                % Medium boundary detection - stepping method
                [x_new, y_new, z_new, crossed_boundary, boundary_info] = movePhotonMediumBoundary(...
                    x, y, z, P, dd_current, lsa, MediaTypeMapping, Size_per_Pixel, X, Y, L, current_tissueType);
                
            case 3
                % Complex boundary detection - multiple detection
                [x_new, y_new, z_new, crossed_boundary, boundary_info] = movePhotonComplexBoundary(...
                    x, y, z, P, dd_current, lsa, MediaTypeMapping, Size_per_Pixel, X, Y, L, current_tissueType);
                
            otherwise
                % Default to simple boundary detection
                [x_new, y_new, z_new, crossed_boundary, boundary_info] = movePhotonSimpleBoundary(...
                    x, y, z, P, dd_current, lsa, MediaTypeMapping, Size_per_Pixel, X, Y, L, current_tissueType);
        end
        
        % Update position
        x = x_new;
        y = y_new;
        z = z_new;
        
        % Calculate actual movement distance
        dx = x - x_prev;
        dy = y - y_prev;
        dz = z - z_prev;
        s_actual = sqrt(dx^2 + dy^2 + dz^2);
        
        % Update time
        t = t + s_actual * (n_water / c0);
        
        % Process g1 and absorption for boundary crossing
        if crossed_boundary
            % Process first segment path (in original medium)
            if boundary_info.segment1_distance > 0
                % Calculate g1 accumulation for first segment
                voxel_x = floor(x_prev / Size_per_Pixel + X/2 + 1);
                voxel_y = floor(y_prev / Size_per_Pixel + Y/2 + 1);
                voxel_z = floor(z_prev / zSize + 1);
                
                if voxel_x >= 1 && voxel_x <= X && voxel_y >= 1 && voxel_y <= Y && voxel_z >= 1 && voxel_z <= L
                    vel_x = VelocityMappingx(voxel_x, voxel_y, voxel_z);
                    vel_y = VelocityMappingy(voxel_x, voxel_y, voxel_z);
                    vel_z = VelocityMappingz(voxel_x, voxel_y, voxel_z);
                    
                    q_segment = P(3,:) * boundary_info.segment1_distance;
                    dphase_segment = dot(q_segment, [vel_x, vel_y, vel_z]);
                    dphase = dphase + dphase_segment;
                    
                    g1 = exp(1j*(2*pi/lambda)*dphase);
                end
                
                % Apply absorption for first segment
                weight = weight * exp(-ua(current_tissueType) * boundary_info.segment1_distance);
            end
            
            % Process second segment path (in new medium)
            if boundary_info.segment2_distance > 0
                segment2_tissue = boundary_info.new_tissueType;
                
                % Update current tissue type
                current_tissueType = segment2_tissue;
                
                % Calculate g1 accumulation for second segment
                mid_x = x_prev + boundary_info.segment1_distance * P(3,1);
                mid_y = y_prev + boundary_info.segment1_distance * P(3,2);
                mid_z = z_prev + boundary_info.segment1_distance * P(3,3);
                
                voxel_x = floor(mid_x / Size_per_Pixel + X/2 + 1);
                voxel_y = floor(mid_y / Size_per_Pixel + Y/2 + 1);
                voxel_z = floor(mid_z / zSize + 1);
                
                if voxel_x >= 1 && voxel_x <= X && voxel_y >= 1 && voxel_y <= Y && voxel_z >= 1 && voxel_z <= L
                    vel_x = VelocityMappingx(voxel_x, voxel_y, voxel_z);
                    vel_y = VelocityMappingy(voxel_x, voxel_y, voxel_z);
                    vel_z = VelocityMappingz(voxel_x, voxel_y, voxel_z);
                    
                    q_segment = P(3,:) * boundary_info.segment2_distance;
                    dphase_segment = dot(q_segment, [vel_x, vel_y, vel_z]);
                    dphase = dphase + dphase_segment;
                    
                    g1 = exp(1j*(2*pi/lambda)*dphase);
                end
                
                % Apply absorption for second segment
                weight = weight * exp(-ua(segment2_tissue) * boundary_info.segment2_distance);
            end
            
        else
            % No boundary crossing, normal processing of entire path
            % Calculate g1 accumulation
            voxel_x = floor(x_prev / Size_per_Pixel + X/2 + 1);
            voxel_y = floor(y_prev / Size_per_Pixel + Y/2 + 1);
            voxel_z = floor(z_prev / zSize + 1);
            
            if voxel_x >= 1 && voxel_x <= X && voxel_y >= 1 && voxel_y <= Y && voxel_z >= 1 && voxel_z <= L
                vel_x = VelocityMappingx(voxel_x, voxel_y, voxel_z);
                vel_y = VelocityMappingy(voxel_x, voxel_y, voxel_z);
                vel_z = VelocityMappingz(voxel_x, voxel_y, voxel_z);
                
                q_total = P(3,:) * s_actual;
                dphase_total = dot(q_total, [vel_x, vel_y, vel_z]);
                dphase = dphase + dphase_total;
                
                g1 = exp(1j*(2*pi/lambda)*dphase);
            end
            
            % Apply absorption
            weight = weight * exp(-ua(current_tissueType) * s_actual);
        end
        
        % ==================== Layer boundary detection and recording ====================
        
        % Reflection detection (upward movement)
        if P(3, 3) < 0 % Upward movement
            if z_prev > z % z decreasing (photon moving upward)
                % Find all integer layers crossed during this movement
                z_start = z_prev/zSize;
                z_end = z/zSize;
                
                % Layers crossed (from high to low)
                layers_crossed = ceil(z_end):floor(z_start);
                layers_crossed = layers_crossed(layers_crossed >= 0 & layers_crossed <= L);
                
                if ~isempty(layers_crossed)
                    for k = layers_crossed
                        % If k is exactly at z_prev (no crossing), skip
                        if k == z_prev
                            continue;
                        end
                        
                        % Calculate intersection with layer boundary
                        fraction = (k - z_prev) / (z - z_prev);
                        x_k = x_prev + fraction * (x - x_prev);
                        y_k = y_prev + fraction * (y - y_prev);
                        z_k = k;
                        
                        % Detection direction for reflection (upward)
                        ud = 0; vd = 0; wd = -1; % Detection direction (toward -z)
                        mu = P(3,1)*ud + P(3,2)*vd + P(3,3)*wd; % Dot product
                        nu = sqrt(1 - mu*mu); % Sine of angle
                        
                        % Calculate angle between photon direction and optical axis (-z direction)
                        photon_direction = P(3, :); % Photon direction vector
                        optical_axis = [0, 0, -1]; % Optical axis direction (-z)
                        cos_angle = dot(photon_direction, optical_axis) / (norm(photon_direction) * norm(optical_axis));
                        
                        % Check if within lens NA range
                        if cos_angle >= cos_max_angle
                            % Get current tissue type (at intersection)
                            [indx_k, indy_k, indz_k] = tissueTypeMapping(x_k/Size_per_Pixel, y_k/Size_per_Pixel, z_k/zSize, X, Y, L);
                            tissue_at_k = MediaTypeMapping(indx_k, indy_k, indz_k);
                            
                            % Calculate reflected field
                            if abs(1 - mu) < 1e-11
                                % Photon already moving in detection direction
                                Q = P;
                                theta = acos(1);
                                pInd = floor(theta/pi*(NN-1)) + 1;
                                s1V = s1tab(tissue_at_k, pInd);
                                s2V = s2tab(tissue_at_k, pInd);
                                s2sq = abs(s2V)^2;
                                s1sq = abs(s1V)^2;
                                F = (s1sq + s2sq)/2;
                                Ed1 = E1 * s2V / sqrt(F);
                                Ed2 = E2 * s1V / sqrt(F);
                                
                            elseif abs(1 + mu) < 1e-11
                                % Photon opposite to detection direction
                                Q = [P(1,:); -P(2,:); -P(3,:)];
                                theta = acos(-1);
                                pInd = floor(theta/pi*(NN-1)) + 1;
                                s1V = s1tab(tissue_at_k, pInd);
                                s2V = s2tab(tissue_at_k, pInd);
                                s2sq = abs(s2V)^2;
                                s1sq = abs(s1V)^2;
                                F = (s1sq + s2sq)/2;
                                Ed1 = E1 * s2V / sqrt(F);
                                Ed2 = E2 * s1V / sqrt(F);
                                
                            else
                                % General case
                                % Calculate rotation axis
                                pVec = [(P(3,2)*wd - P(3,3)*vd), ...
                                    (P(3,3)*ud - P(3,1)*wd), ...
                                    (P(3,1)*vd - P(3,2)*ud)] / nu;
                                
                                % Calculate rotation angle
                                cosphi = P(2,1)*pVec(1) + P(2,2)*pVec(2) + P(2,3)*pVec(3);
                                sinphi = -(P(1,1)*pVec(1) + P(1,2)*pVec(2) + P(1,3)*pVec(3));
                                
                                % Rotation matrix
                                A = [mu*cosphi, mu*sinphi, -nu; ...
                                    -sinphi,   cosphi,    0; ...
                                    nu*cosphi, nu*sinphi, mu];
                                
                                % Apply rotation
                                Q = A * P;
                                
                                % Get scattering amplitudes
                                theta = acos(mu);
                                pInd = floor(theta/pi*(NN-1)) + 1;
                                s1V = s1tab(tissue_at_k, pInd);
                                s2V = s2tab(tissue_at_k, pInd);
                                s2sq = abs(s2V)^2;
                                s1sq = abs(s1V)^2;
                                
                                % Calculate field components
                                e1sq = abs(E1)^2;
                                e2sq = abs(E2)^2;
                                e12 = real(E1*conj(E2));
                                F = (s2sq*e1sq + s1sq*e2sq)*cosphi^2 + ...
                                    (s1sq*e1sq + s2sq*e2sq)*sinphi^2 + ...
                                    2*(s2sq - s1sq)*e12*cosphi*sinphi;
                                
                                Ed1 = (cosphi*E1 + sinphi*E2)*s2V/sqrt(F);
                                Ed2 = (-sinphi*E1 + cosphi*E2)*s1V/sqrt(F);
                            end
                            
                            % Calculate electric field components
                            Ex = Ed1*Q(1,1) + Ed2*Q(2,1);
                            Ey = Ed1*Q(1,2) + Ed2*Q(2,2);
                            % Calculate total electric field magnitude (non-polarized)
                            E_field = sqrt(abs(Ex)^2 + abs(Ey)^2) * exp(1j*(angle(Ex) + angle(Ey)));
                            
                            % Determine layer index (k is boundary, layer = k+1)
                            layer_index = floor(z_k) + 1;
                            if layer_index < 1 || layer_index > L
                                continue; % Skip invalid layer
                            end
                            
                            % Calculate pixel coordinates
                            RecordX = round( (x_k/Size_per_Pixel + X/2) );
                            RecordY = round( (y_k/Size_per_Pixel + Y/2) );
                            
                            % Ensure pixel coordinates are within range
                            if RecordX >= 1 && RecordX <= scaled_X && RecordY >= 1 && RecordY <= scaled_Y
                                % Get current g1 value (considering time scale)
                                g1_value = g1^(timeScale);
                                
                                % Process for each polarization angle
                                for pol_idx = 1:num_pol_angles
                                    pol_angle = pol_angles(pol_idx);
                                    
                                    % Calculate polarization direction vector
                                    pol_vector = [cosd(pol_angle); sind(pol_angle); 0];
                                    
                                    % Calculate electric field projection onto polarization direction
                                    E_proj = Ex * pol_vector(1) + Ey * pol_vector(2);
                                    
                                    % Determine scattering type and record
                                    if nsct - 1 == 1 % Single scattering
                                        IMRS_layer_pol(RecordX, RecordY, layer_index, pol_idx) = ...
                                            IMRS_layer_pol(RecordX, RecordY, layer_index, pol_idx) + E_proj * sqrt(weight);
                                        G1_perm_RS_layer(RecordX, RecordY, layer_index) = ...
                                            G1_perm_RS_layer(RecordX, RecordY, layer_index) + weight * g1_value;
                                        Photonpath_perm_RS_layer(RecordX, RecordY, layer_index) = ...
                                            Photonpath_perm_RS_layer(RecordX, RecordY, layer_index) + 1;
                                    else % Multiple scattering
                                        IMRM_layer_pol(RecordX, RecordY, layer_index, pol_idx) = ...
                                            IMRM_layer_pol(RecordX, RecordY, layer_index, pol_idx) + E_proj * sqrt(weight);
                                        G1_perm_RM_layer(RecordX, RecordY, layer_index) = ...
                                            G1_perm_RM_layer(RecordX, RecordY, layer_index) + weight * g1_value;
                                        Photonpath_perm_RM_layer(RecordX, RecordY, layer_index) = ...
                                            Photonpath_perm_RM_layer(RecordX, RecordY, layer_index) + 1;
                                    end
                                end
                                
                                % Record non-polarized image
                                if nsct - 1 == 1
                                    IMRS_layer(RecordX, RecordY, layer_index) = ...
                                        IMRS_layer(RecordX, RecordY, layer_index) + E_field * sqrt(weight);
                                else
                                    IMRM_layer(RecordX, RecordY, layer_index) = ...
                                        IMRM_layer(RecordX, RecordY, layer_index) + E_field * sqrt(weight);
                                end
                            end
                        end
                    end
                end
            end
        end
        
        % Transmission detection (downward movement)
        if P(3, 3) > 0 % Downward movement
            if z_prev < z % z increasing (photon moving downward)
                % Find all integer layers crossed during this movement
                z_start = z_prev/zSize;
                z_end = z/zSize;
                
                % Layers crossed (from low to high)
                layers_crossed = ceil(z_start):floor(z_end);
                layers_crossed = layers_crossed(layers_crossed >= 0 & layers_crossed <= L);
                
                if ~isempty(layers_crossed)
                    for k = layers_crossed
                        if k == z_prev
                            continue;
                        end
                        
                        % Calculate intersection with layer boundary
                        fraction = (k - z_prev) / (z - z_prev);
                        x_k = x_prev + fraction * (x - x_prev);
                        y_k = y_prev + fraction * (y - y_prev);
                        z_k = k;
                        
                        % Detection direction for transmission (downward)
                        ud = 0; vd = 0; wd = 1; % Detection direction (toward +z)
                        mu = P(3,1)*ud + P(3,2)*vd + P(3,3)*wd;
                        nu = sqrt(1 - mu*mu);
                        
                        % Calculate angle between photon direction and optical axis (+z direction)
                        photon_direction = P(3, :);
                        optical_axis = [0, 0, 1];
                        cos_angle = dot(photon_direction, optical_axis) / (norm(photon_direction) * norm(optical_axis));
                        
                        % Check if within lens NA range
                        if cos_angle >= cos_max_angle
                            % Get current tissue type (at intersection)
                            [indx_k, indy_k, indz_k] = tissueTypeMapping(x_k/Size_per_Pixel, y_k/Size_per_Pixel, z_k/zSize, X, Y, L);
                            tissue_at_k = MediaTypeMapping(indx_k, indy_k, indz_k);
                            
                            % Calculate transmitted field (similar to reflection calculation)
                            if abs(1 - mu) < 1e-11
                                Q = P;
                                theta = acos(1);
                                pInd = floor(theta/pi*(NN-1)) + 1;
                                s1V = s1tab(tissue_at_k, pInd);
                                s2V = s2tab(tissue_at_k, pInd);
                                s2sq = abs(s2V)^2;
                                s1sq = abs(s1V)^2;
                                F = (s1sq + s2sq)/2;
                                Ed1 = E1 * s2V / sqrt(F);
                                Ed2 = E2 * s1V / sqrt(F);
                                
                            elseif abs(1 + mu) < 1e-11
                                Q = [P(1,:); -P(2,:); -P(3,:)];
                                theta = acos(-1);
                                pInd = floor(theta/pi*(NN-1)) + 1;
                                s1V = s1tab(tissue_at_k, pInd);
                                s2V = s2tab(tissue_at_k, pInd);
                                s2sq = abs(s2V)^2;
                                s1sq = abs(s1V)^2;
                                F = (s1sq + s2sq)/2;
                                Ed1 = E1 * s2V / sqrt(F);
                                Ed2 = E2 * s1V / sqrt(F);
                                
                            else
                                pVec = [(P(3,2)*wd - P(3,3)*vd), ...
                                    (P(3,3)*ud - P(3,1)*wd), ...
                                    (P(3,1)*vd - P(3,2)*ud)] / nu;
                                
                                cosphi = P(2,1)*pVec(1) + P(2,2)*pVec(2) + P(2,3)*pVec(3);
                                sinphi = -(P(1,1)*pVec(1) + P(1,2)*pVec(2) + P(1,3)*pVec(3));
                                
                                A = [mu*cosphi, mu*sinphi, -nu; ...
                                    -sinphi,   cosphi,    0; ...
                                    nu*cosphi, nu*sinphi, mu];
                                
                                Q = A * P;
                                
                                theta = acos(mu);
                                pInd = floor(theta/pi*(NN-1)) + 1;
                                s1V = s1tab(tissue_at_k, pInd);
                                s2V = s2tab(tissue_at_k, pInd);
                                s2sq = abs(s2V)^2;
                                s1sq = abs(s1V)^2;
                                
                                e1sq = abs(E1)^2;
                                e2sq = abs(E2)^2;
                                e12 = real(E1*conj(E2));
                                F = (s2sq*e1sq + s1sq*e2sq)*cosphi^2 + ...
                                    (s1sq*e1sq + s2sq*e2sq)*sinphi^2 + ...
                                    2*(s2sq - s1sq)*e12*cosphi*sinphi;
                                
                                Ed1 = (cosphi*E1 + sinphi*E2)*s2V/sqrt(F);
                                Ed2 = (-sinphi*E1 + cosphi*E2)*s1V/sqrt(F);
                            end
                            
                            Ex = Ed1*Q(1,1) + Ed2*Q(2,1);
                            Ey = Ed1*Q(1,2) + Ed2*Q(2,2);
                            % Calculate total electric field magnitude (non-polarized)
                            E_field = sqrt(abs(Ex)^2 + abs(Ey)^2) * exp(1j*(angle(Ex) + angle(Ey)));
                            
                            layer_index = floor(z_k) + 1;
                            if layer_index < 1 || layer_index > L
                                continue;
                            end
                            
                            RecordX = round( (x_k/Size_per_Pixel + X/2) );
                            RecordY = round( (y_k/Size_per_Pixel + Y/2) );
                            
                            if RecordX >= 1 && RecordX <= scaled_X && RecordY >= 1 && RecordY <= scaled_Y
                                g1_value = g1^(timeScale);
                                
                                for pol_idx = 1:num_pol_angles
                                    pol_angle = pol_angles(pol_idx);
                                    pol_vector = [cosd(pol_angle); sind(pol_angle); 0];
                                    E_proj = Ex * pol_vector(1) + Ey * pol_vector(2);
                                    
                                    if nsct - 1 == 1
                                        IMTS_layer_pol(RecordX, RecordY, layer_index, pol_idx) = ...
                                            IMTS_layer_pol(RecordX, RecordY, layer_index, pol_idx) + E_proj * sqrt(weight);
                                        G1_perm_TS_layer(RecordX, RecordY, layer_index) = ...
                                            G1_perm_TS_layer(RecordX, RecordY, layer_index) + weight * g1_value;
                                        Photonpath_perm_TS_layer(RecordX, RecordY, layer_index) = ...
                                            Photonpath_perm_TS_layer(RecordX, RecordY, layer_index) + 1;
                                    else
                                        IMTM_layer_pol(RecordX, RecordY, layer_index, pol_idx) = ...
                                            IMTM_layer_pol(RecordX, RecordY, layer_index, pol_idx) + E_proj * sqrt(weight);
                                        G1_perm_TM_layer(RecordX, RecordY, layer_index) = ...
                                            G1_perm_TM_layer(RecordX, RecordY, layer_index) + weight * g1_value;
                                        Photonpath_perm_TM_layer(RecordX, RecordY, layer_index) = ...
                                            Photonpath_perm_TM_layer(RecordX, RecordY, layer_index) + 1;
                                    end
                                end
                                
                                % Record non-polarized image
                                if nsct - 1 == 1
                                    IMTS_layer(RecordX, RecordY, layer_index) = ...
                                        IMTS_layer(RecordX, RecordY, layer_index) + E_field * sqrt(weight);
                                else
                                    IMTM_layer(RecordX, RecordY, layer_index) = ...
                                        IMTM_layer(RecordX, RecordY, layer_index) + E_field * sqrt(weight);
                                end
                            end
                        end
                    end
                end
            end
        end
        
        % ==================== Surface detection (for termination) ====================
        
        % Check if photon exits medium - transmission (bottom)
        if z >= L * zSize
            % Calculate angle between photon direction and optical axis (+z direction)
            photon_direction = P(3, :);
            optical_axis = [0, 0, 1];
            cos_angle = dot(photon_direction, optical_axis) / (norm(photon_direction) * norm(optical_axis));
            
            if cos_angle >= cos_max_angle
                % Photon exits from bottom (transmission)
                break;
            end
            
        % Check if photon exits medium - reflection (top)
        elseif z <= 0
            photon_direction = P(3, :);
            optical_axis = [0, 0, -1];
            cos_angle = dot(photon_direction, optical_axis) / (norm(photon_direction) * norm(optical_axis));
            
            if cos_angle >= cos_max_angle
                % Photon exits from top (reflection)
                break;
            end
        end
        
        % Check photon termination condition
        if weight <= 0.001
            break;
        end
        
        % ==================== Scattering event ====================
        e1sq = abs(E1)^2;
        e2sq = abs(E2)^2;
        e12 = real(E1*conj(E2));
        
        % Rejection sampling for scattering angle mu
        while true
            mu = rand(1)*2 - 1;
            theta = acos(mu);
            pInd = floor(theta/pi*(NN-1)) + 1;
            s1 = s1tab(current_tissueType, pInd);
            s2 = s2tab(current_tissueType, pInd);
            s2sq = abs(s2)^2;
            s1sq = abs(s1)^2;
            rejmu = rejmu + 1;
            
            if rand(1)*psca0(current_tissueType) < (s1sq + s2sq)/2
                break;
            end
        end
        
        % Calculate maximum possible F value
        a = s2sq*e1sq + s1sq*e2sq;
        b = s1sq*e1sq + s2sq*e2sq;
        c = 2*(s2sq - s1sq)*e12;
        Fmax = (a + b)/2 + sqrt((a - b)^2 + c^2)/2;
        
        % Rejection sampling for azimuthal angle phi
        while true
            phi = rand(1)*2*pi;
            cosphi = cos(phi);
            sinphi = sin(phi);
            F = a*cosphi^2 + b*sinphi^2 + c*cosphi*sinphi;
            rejphi = rejphi + 1;
            
            if rand(1)*Fmax < F
                break;
            end
        end
        
        % Update direction after scattering
        nu = sqrt(1 - mu^2);
        A = [mu*cosphi, mu*sinphi, -nu; ...
            -sinphi,   cosphi,    0; ...
            nu*cosphi, nu*sinphi, mu];
        P = A * P;
        
        % Update polarization after scattering
        F = sqrt(F);
        e = (cosphi*E1 + sinphi*E2)*s2/F;
        E2 = (-sinphi*E1 + cosphi*E2)*s1/F;
        E1 = e;
        
        % Russian roulette for photon termination
        if weight < 0.001
            if rand(1) > 0.1
                weight = 0;
            else
                weight = weight / 0.1;
            end
        end
        
        nsct = nsct + 1;
    end
end

% ==================== 4. Store output results ====================

% Store the 16 output variables in iImage structure
iImage.IMRS_layer_pol = IMRS_layer_pol;
iImage.IMRM_layer_pol = IMRM_layer_pol;
iImage.IMRS_layer = IMRS_layer;
iImage.IMRM_layer = IMRM_layer;
iImage.G1_perm_RS_layer = G1_perm_RS_layer;
iImage.G1_perm_RM_layer = G1_perm_RM_layer;
iImage.Photonpath_perm_RS_layer = Photonpath_perm_RS_layer;
iImage.Photonpath_perm_RM_layer = Photonpath_perm_RM_layer;

iImage.IMTS_layer_pol = IMTS_layer_pol;
iImage.IMTM_layer_pol = IMTM_layer_pol;
iImage.IMTS_layer = IMTS_layer;
iImage.IMTM_layer = IMTM_layer;
iImage.G1_perm_TS_layer = G1_perm_TS_layer;
iImage.G1_perm_TM_layer = G1_perm_TM_layer;
iImage.Photonpath_perm_TS_layer = Photonpath_perm_TS_layer;
iImage.Photonpath_perm_TM_layer = Photonpath_perm_TM_layer;

% Store additional information
iImage.scaled_X = scaled_X;
iImage.scaled_Y = scaled_Y;
iImage.pol_angles = pol_angles;

toc;
disp('Photon transport completed');
end

% ==================== Boundary detection helper functions ====================

function [x_new, y_new, z_new, crossed_boundary, boundary_info] = movePhotonNoBoundaryCheck(...
    x, y, z, P, dd, lsa_current, current_tissueType)
% Homogeneous media, no boundary detection
x_new = x + dd * P(3, 1) * lsa_current;
y_new = y + dd * P(3, 2) * lsa_current;
z_new = z + dd * P(3, 3) * lsa_current;
crossed_boundary = false;
boundary_info = struct();
end

function [x_new, y_new, z_new, crossed_boundary, boundary_info] = movePhotonSimpleBoundary(...
    x, y, z, P, dd, lsa, MediaTypeMapping, Size_per_Pixel, X, Y, L, current_tissueType)
% Simple boundary detection - bisection method

x_prev = x; y_prev = y; z_prev = z;
total_distance = dd * lsa(current_tissueType);

% Calculate target position
x_target = x + total_distance * P(3, 1);
y_target = y + total_distance * P(3, 2);
z_target = z + total_distance * P(3, 3);

% Check tissue type at target position
[indx_target, indy_target, indz_target] = tissueTypeMapping(...
    x_target/Size_per_Pixel, y_target/Size_per_Pixel, z_target/Size_per_Pixel, X, Y, L);
target_tissueType = MediaTypeMapping(indx_target, indy_target, indz_target);

crossed_boundary = (target_tissueType ~= current_tissueType);
boundary_info = struct();

if ~crossed_boundary
    % No boundary crossing, move directly to target position
    x_new = x_target;
    y_new = y_target;
    z_new = z_target;
else
    % Use bisection method to find boundary point
    low = 0;
    high = total_distance;
    tolerance = Size_per_Pixel * 0.01;
    
    boundary_found = false;
    boundary_distance = total_distance;
    
    for iter = 1:20
        mid = (low + high) / 2;
        x_mid = x + mid * P(3, 1);
        y_mid = y + mid * P(3, 2);
        z_mid = z + mid * P(3, 3);
        
        [indx_mid, indy_mid, indz_mid] = tissueTypeMapping(...
            x_mid/Size_per_Pixel, y_mid/Size_per_Pixel, z_mid/Size_per_Pixel, X, Y, L);
        mid_tissueType = MediaTypeMapping(indx_mid, indy_mid, indz_mid);
        
        if mid_tissueType == current_tissueType
            low = mid;
        else
            high = mid;
            boundary_distance = mid;
            boundary_found = true;
        end
        
        if (high - low) < tolerance
            break;
        end
    end
    
    if boundary_found
        % Move to boundary point
        x_new = x + boundary_distance * P(3, 1);
        y_new = y + boundary_distance * P(3, 2);
        z_new = z + boundary_distance * P(3, 3);
        
        % Get tissue type of new medium
        [indx_new, indy_new, indz_new] = tissueTypeMapping(...
            x_new/Size_per_Pixel, y_new/Size_per_Pixel, z_new/Size_per_Pixel, X, Y, L);
        new_tissueType = MediaTypeMapping(indx_new, indy_new, indz_new);
        
        boundary_info.segment1_distance = boundary_distance;
        boundary_info.segment2_distance = total_distance - boundary_distance;
        boundary_info.new_tissueType = new_tissueType;
    else
        % Boundary not found, use target position
        x_new = x_target;
        y_new = y_target;
        z_new = z_target;
        crossed_boundary = false;
    end
end
end

function [x_new, y_new, z_new, crossed_boundary, boundary_info] = movePhotonMediumBoundary(...
    x, y, z, P, dd, lsa, MediaTypeMapping, Size_per_Pixel, X, Y, L, current_tissueType)
% Medium boundary detection - stepping method

x_prev = x; y_prev = y; z_prev = z;
total_distance = dd * lsa(current_tissueType);
step_size = Size_per_Pixel * 0.5;

current_pos = [x, y, z];
direction = P(3, :);
remaining_distance = total_distance;
current_tissue = current_tissueType;

boundary_crossed = false;
boundary_distance = 0;

% Step movement
while remaining_distance > 0 && ~boundary_crossed
    step = min(step_size, remaining_distance);
    next_pos = current_pos + step * direction;
    
    % Check tissue type
    [indx_next, indy_next, indz_next] = tissueTypeMapping(...
        next_pos(1)/Size_per_Pixel, next_pos(2)/Size_per_Pixel, next_pos(3)/Size_per_Pixel, X, Y, L);
    next_tissue = MediaTypeMapping(indx_next, indy_next, indz_next);
    
    if next_tissue == current_tissue
        current_pos = next_pos;
        remaining_distance = remaining_distance - step;
        boundary_distance = boundary_distance + step;
    else
        boundary_crossed = true;
        
        % Use bisection to precisely locate boundary point
        low_pos = current_pos;
        high_pos = next_pos;
        
        for refine = 1:10
            mid_pos = (low_pos + high_pos) / 2;
            [indx_mid, indy_mid, indz_mid] = tissueTypeMapping(...
                mid_pos(1)/Size_per_Pixel, mid_pos(2)/Size_per_Pixel, mid_pos(3)/Size_per_Pixel, X, Y, L);
            mid_tissue = MediaTypeMapping(indx_mid, indy_mid, indz_mid);
            
            if mid_tissue == current_tissue
                low_pos = mid_pos;
            else
                high_pos = mid_pos;
            end
        end
        
        current_pos = low_pos;
        boundary_distance = norm(current_pos - [x, y, z]);
    end
end

x_new = current_pos(1);
y_new = current_pos(2);
z_new = current_pos(3);

boundary_info = struct();
if boundary_crossed
    crossed_boundary = true;
    [indx_new, indy_new, indz_new] = tissueTypeMapping(...
        x_new/Size_per_Pixel, y_new/Size_per_Pixel, z_new/Size_per_Pixel, X, Y, L);
    new_tissueType = MediaTypeMapping(indx_new, indy_new, indz_new);
    
    boundary_info.segment1_distance = boundary_distance;
    boundary_info.segment2_distance = total_distance - boundary_distance;
    boundary_info.new_tissueType = new_tissueType;
else
    crossed_boundary = false;
end
end

function [x_new, y_new, z_new, crossed_boundary, boundary_info] = movePhotonComplexBoundary(...
    x, y, z, P, dd, lsa, MediaTypeMapping, Size_per_Pixel, X, Y, L, current_tissueType)
% Complex boundary detection - multiple detection and step scattering

x_prev = x; y_prev = y; z_prev = z;
total_distance = dd * lsa(current_tissueType);
step_size = Size_per_Pixel * 0.1;

current_pos = [x, y, z];
direction = P(3, :);
remaining_distance = total_distance;
current_tissue = current_tissueType;

segments = [];
current_segment_distance = 0;
current_segment_tissue = current_tissue;

% Detailed step detection
while remaining_distance > 0
    step = min(step_size, remaining_distance);
    next_pos = current_pos + step * direction;
    
    % Check tissue type
    [indx_next, indy_next, indz_next] = tissueTypeMapping(...
        next_pos(1)/Size_per_Pixel, next_pos(2)/Size_per_Pixel, next_pos(3)/Size_per_Pixel, X, Y, L);
    next_tissue = MediaTypeMapping(indx_next, indy_next, indz_next);
    
    if next_tissue == current_segment_tissue
        current_segment_distance = current_segment_distance + step;
    else
        if current_segment_distance > 0
            segments(end+1) = struct(...
                'distance', current_segment_distance, ...
                'tissueType', current_segment_tissue);
        end
        
        current_segment_tissue = next_tissue;
        current_segment_distance = step;
    end
    
    current_pos = next_pos;
    remaining_distance = remaining_distance - step;
end

% Record last segment
if current_segment_distance > 0
    segments(end+1) = struct(...
        'distance', current_segment_distance, ...
        'tissueType', current_segment_tissue);
end

x_new = current_pos(1);
y_new = current_pos(2);
z_new = current_pos(3);

if length(segments) > 1
    crossed_boundary = true;
    boundary_info.segments = segments;
    boundary_info.num_segments = length(segments);
    boundary_info.segment1_distance = segments(1).distance;
    boundary_info.segment2_distance = total_distance - segments(1).distance;
    boundary_info.new_tissueType = segments(2).tissueType;
else
    crossed_boundary = false;
    boundary_info = struct();
end
end