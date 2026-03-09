%%

Modelsetting;

%%
num_runs = 10;
num_core = 10;

% Check if parallel pool is open, if not, open it
if isempty(gcp('nocreate'))
    parpool(num_core);
end

% Preallocate result storage
all_iImages = cell(num_runs, 1);

fprintf('Starting %d parallel simulations...\n', num_runs);

% Execute multiple simulations in parallel
parfor run_idx = 1:num_runs
    fprintf('Process %d: Starting simulation %d\n', labindex, run_idx);
    
    % Create independent random number stream for each run to ensure different results
    rng_stream = RandStream('mlfg6331_64', 'Seed', run_idx);
    savedStream = RandStream.getGlobalStream();
    RandStream.setGlobalStream(rng_stream);
    
    try
        % Execute simulation
        [MieC, iImage_temp] = dEMC(MieScatter, Light, Dynamic, G, Camera);
        
        % Store results
        all_iImages{run_idx} = iImage_temp;
        
        fprintf('Process %d: Simulation %d completed\n', labindex, run_idx);
        
    catch ME
        fprintf('Process %d: Simulation %d failed: %s\n', labindex, run_idx, ME.message);
        % If failed, use empty results
        all_iImages{run_idx} = [];
    end
    
    % Restore original random number stream
    RandStream.setGlobalStream(savedStream);
end

% Restore serial random number stream
rng('shuffle');

% Filter out failed results
valid_indices = ~cellfun(@isempty, all_iImages);
all_iImages = all_iImages(valid_indices);

if isempty(all_iImages)
    error('All simulation runs failed!');
end

fprintf('Successfully completed %d simulations, starting to combine results...\n', length(all_iImages));


% Initialize combined result structure
first_iImage = all_iImages{1};
Image = struct();

% List of all image fields to combine (including G1 and Photonpath fields)
image_fields = {'IMRS_layer_pol', 'IMRM_layer_pol', 'IMRS_layer', 'IMRM_layer', ...
               'G1_perm_RS_layer', 'G1_perm_RM_layer', 'Photonpath_perm_RS_layer', 'Photonpath_perm_RM_layer', ...
               'IMTS_layer_pol', 'IMTM_layer_pol', 'IMTS_layer', 'IMTM_layer', ...
               'G1_perm_TS_layer', 'G1_perm_TM_layer', 'Photonpath_perm_TS_layer', 'Photonpath_perm_TM_layer'};

% Copy non-image fields
non_image_fields = {'scaled_X', 'scaled_Y', 'pol_angles'};
for i = 1:length(non_image_fields)
    field = non_image_fields{i};
    if isfield(first_iImage, field)
        Image.(field) = first_iImage.(field);
    end
end

% Combine image fields by summation
for i = 1:length(image_fields)
    field_name = image_fields{i};
    if isfield(first_iImage, field_name)
        fprintf('Combining field: %s\n', field_name);
        
        % Get dimension information of the field
        field_data = first_iImage.(field_name);
        dims = size(field_data);
        ndims_field = ndims(field_data);
        
        % Initialize sum result
        if isreal(field_data)
            sum_field = zeros(dims);
        else
            sum_field = complex(zeros(dims));
        end
        
        % Sum results from all runs
        for j = 1:length(all_iImages)
            if ~isempty(all_iImages{j}) && isfield(all_iImages{j}, field_name)
                sum_field = sum_field + all_iImages{j}.(field_name);
            end
        end
        
        % Store combined result
        Image.(field_name) = sum_field;
    end
end

% Calculate normalized G1 values by dividing G1_perm by Photonpath_perm
fprintf('Calculating normalized G1 values...\n');

% Reflection: Single scattering
if isfield(Image, 'G1_perm_RS_layer') && isfield(Image, 'Photonpath_perm_RS_layer')
    valid_mask = Image.Photonpath_perm_RS_layer > 0;
    Image.G1RS_layer = zeros(size(Image.G1_perm_RS_layer));
    Image.G1RS_layer(valid_mask) = Image.G1_perm_RS_layer(valid_mask) ./ ...
                                         Image.Photonpath_perm_RS_layer(valid_mask);
    fprintf('Calculated G1RS_layer\n');
end

% Reflection: Multiple scattering
if isfield(Image, 'G1_perm_RM_layer') && isfield(Image, 'Photonpath_perm_RM_layer')
    valid_mask = Image.Photonpath_perm_RM_layer > 0;
    Image.G1RM_layer = zeros(size(Image.G1_perm_RM_layer));
    Image.G1RM_layer(valid_mask) = Image.G1_perm_RM_layer(valid_mask) ./ ...
                                         Image.Photonpath_perm_RM_layer(valid_mask);
    fprintf('Calculated G1RM_layer\n');
end

% Transmission: Single scattering
if isfield(Image, 'G1_perm_TS_layer') && isfield(Image, 'Photonpath_perm_TS_layer')
    valid_mask = Image.Photonpath_perm_TS_layer > 0;
    Image.G1TS_layer = zeros(size(Image.G1_perm_TS_layer));
    Image.G1TS_layer(valid_mask) = Image.G1_perm_TS_layer(valid_mask) ./ ...
                                         Image.Photonpath_perm_TS_layer(valid_mask);
    fprintf('Calculated G1TS_layer\n');
end

% Transmission: Multiple scattering
if isfield(Image, 'G1_perm_TM_layer') && isfield(Image, 'Photonpath_perm_TM_layer')
    valid_mask = Image.Photonpath_perm_TM_layer > 0;
    Image.G1TM_layer = zeros(size(Image.G1_perm_TM_layer));
    Image.G1TM_layer(valid_mask) = Image.G1_perm_TM_layer(valid_mask) ./ ...
                                         Image.Photonpath_perm_TM_layer(valid_mask);
    fprintf('Calculated G1TM_layer\n');
end

%%
RIs = OUprocess (Image.IMRS_layer, Image.G1RS_layer, Camera, Light); % Rflection single scattering
RIm = OUprocess (Image.IMRM_layer, Image.G1RM_layer, Camera, Light); % Rflection multiple scattering
TIs = OUprocess (Image.IMTS_layer, Image.G1TS_layer, Camera, Light); % Transmission single scattering
TIm = OUprocess (Image.IMTM_layer, Image.G1TM_layer, Camera, Light); % Transmission multiple scattering
%% OU process and speckle size modulation
function [Id] = OUprocess (E0, g, Camera, Light)
sig = mean(mean(abs(E0).^2, 1), 2);
snap = Camera.snap;
T = Camera.T;
fps = Camera.fps;
gap = round(1*10^6/fps - T);
[X, Y, Z] = size(E0);
num_time_frames = round(T/snap);

% PreCalculation
Ss = 2.44 * Light.lambda * (1 + Camera.M) * Camera.f_num;
sp = Ss./Camera.xSize;
sigma = sp / (2 * sqrt(2));
kernel_radius = ceil(3 * sigma);
kernel_size = 2 * kernel_radius + 1;
[x, y] = meshgrid(-kernel_radius:kernel_radius, -kernel_radius:kernel_radius);
gaussian_kernel = exp(-(x.^2 + y.^2) / (2 * sigma^2));
gaussian_kernel = gaussian_kernel / sum(gaussian_kernel(:));


g_real = real(g);
g_imag = imag(g);

Id = zeros(X, Y, Z, Nframe, 'like', real(E0));
E_current = E0;

for iter = 1:Nframe
    fprintf('Iteration %d/%d\n', iter, Nframe);
    
   
    E_time_series = zeros(X, Y, Z, num_time_frames, 'like', E0);
    E_time_series(:, :, :, 1) = E_current;
    
    % OU 
    for tt = 2:num_time_frames
        E_prev = E_time_series(:, :, :, tt-1);
        
        sigma2 = sig;
        x1 = real(E_prev);
        y1 = imag(E_prev);
        
        cond_var = max((sigma2 / 2) .* (1 - abs(g).^2), 1e-10);
        mu_x2 = g_real .* x1 - g_imag .* y1;
        mu_y2 = g_imag .* x1 + g_real .* y1;
        
        noise_x = sqrt(cond_var) .* randn(size(E_prev));
        noise_y = sqrt(cond_var) .* randn(size(E_prev));
        
        E_time_series(:, :, :, tt) = (mu_x2 + noise_x) + 1i * (mu_y2 + noise_y);
    end
    
    % speckle size
    for tt = 1:num_time_frames
        E_frame = E_time_series(:, :, :, tt);
        E_real = real(E_frame);
        E_imag = imag(E_frame);
        
        
        for z = 1:Z
            E_real_filtered = imfilter(E_real(:, :, z), gaussian_kernel, 'same', 'conv');
            E_imag_filtered = imfilter(E_imag(:, :, z), gaussian_kernel, 'same', 'conv');
            E_time_series(:, :, z, tt) = E_real_filtered + 1i * E_imag_filtered;
        end
    end
    
   
    I_ca = sum(abs(E_time_series).^2, 4);
    Id(:, :, :, iter) = I_ca;
    
   
    E_current = E_time_series(:, :, :, end);
    if gap > 0
        for tt2 = 1:gap
            % ÄÚÁªElectricFiledCopula¼ÆËã
            sigma2 = sig;
            x1 = real(E_current);
            y1 = imag(E_current);
            
            cond_var = max((sigma2 / 2) .* (1 - abs(g).^2), 1e-10);
            mu_x2 = g_real .* x1 - g_imag .* y1;
            mu_y2 = g_imag .* x1 + g_real .* y1;
            
            noise_x = sqrt(cond_var) .* randn(size(E_current));
            noise_y = sqrt(cond_var) .* randn(size(E_current));
            
            E_current = (mu_x2 + noise_x) + 1i * (mu_y2 + noise_y);
        end
    end
end
end
