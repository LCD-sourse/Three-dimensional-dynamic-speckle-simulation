%% Parameter Setting
% Please set your model parameters before starting

% Author: LI Chengdong in Matrix Optics Research Group
% Shanghai Jiao Tong University

%% MieScatter Parameters
MieScatter = struct();

MieScatter.n_water = 1.33; % Refractive index of water (determined by wavelength)
MieScatter.NN = 10001; % Length of S1 and S2 tables
MieScatter.MieTypeNum = 2; % Number of Mie scatterer types
MieScatter.n_sctCOMPLEX = [1.4 + 0.0000*1i 1.4 + 0.0000*1i]; % Complex refractive index of scatterers (determined by wavelength)
MieScatter.ua = [0*10^(-3) 0*10^(-3)]; % Absorption coefficient (¦Ìm??)
MieScatter.dia = [2 2]; % Diameter of scatterers (¦Ìm)
MieScatter.Volfrac = [0.45 0.45]; % Volume fraction of Mie scatterers

%% Light Source Parameters
Light = struct();
Light.PhotonNum = 1e5; % Number of photons in a single photon packet
Light.c0 = 299.8; % Speed of light in vacuum (¦Ìm/ps)
Light.lambda = 0.785; % Wavelength in vacuum (¦Ìm) ********************  
Light.Xrange = 128*30; % Illumination area in x-direction (¦Ìm)
Light.Yrange = 128*30; % Illumination area in y-direction (¦Ìm)

%% Geometry Parameters
X = 128; Y = 128; L = 25; % Number of voxels in each dimension
G = struct();
% Set media types
G.MediaTypeMapping = ones(X, Y, L); % Initialize with homogeneous medium

% Load vessel data and create media type mapping
load('D:\Data\Vessels.mat');
G.MediaTypeMapping = flip(round(imresize3(Vessel{8}(:,:,end:-1:1), [128,128,25])) + 1, 3);

G.X = X; % Number of voxels in x-direction
G.Y = Y; % Number of voxels in y-direction
G.L = L; % Number of voxels in z-direction
G.complex = 1; % Media complexity level: 
               % 0: Homogeneous medium
               % 1~3: Increasing media complexity. Higher values use more 
               %      sophisticated boundary detection algorithms but 
               %      increase computation time.

%% Dynamic Parameters
Dynamic = struct();

% Brownian motion diffusion coefficient mapping
Dynamic.BMMapping = ones(X, Y, L) * 2.45 * 10^(-13); % m?/s

% Generate random flow directions
theta = rand(X*Y*L, 1) * pi; 
phi = rand(X*Y*L, 1) * 2 * pi;

% Identify dynamic regions (vessel regions)
vol = uint8(reshape(G.MediaTypeMapping, X*Y*L, 1));
v_dynamic = vol == 2; % Assuming type 2 corresponds to vessels

% Initialize velocity components (static regions)
vz = 0 * cos(theta);   % mm/s - ¦Ìm/ms
vx = 0 * sin(theta) .* cos(phi); 
vy = 0 * sin(theta) .* sin(phi); 

% Assign velocity to dynamic regions (vessels)
vz(v_dynamic) = 10 * cos(theta(v_dynamic));   % mm/s - ¦Ìm/ms
vx(v_dynamic) = 10 * sin(theta(v_dynamic)) .* cos(phi(v_dynamic)); 
vy(v_dynamic) = 10 * sin(theta(v_dynamic)) .* sin(phi(v_dynamic)); 

% Reshape to 3D mapping
Dynamic.FVxMapping = reshape(vx, [X, Y, L]);
Dynamic.FVyMapping = reshape(vy, [X, Y, L]);
Dynamic.FVzMapping = reshape(vz, [X, Y, L]);

%% Reflection Camera Parameters
Camera = struct();
Camera.snap = 100; % Sampling/snapshot time (¦Ìs)
Camera.T = 5000;   % Exposure time (¦Ìs)
Camera.fps = 40;   % Frame rate (Hz)
Camera.NA = 1;     % Numerical aperture
Camera.xSize = 40; % Pixel size in x and y directions (¦Ìm)
Camera.zSize = 7;  % Layer thickness in z-direction (¦Ìm)
Camera.pol_angles = [0, 45, 90, 135]; % Polarization angles (degrees)
Camera.M = 10;     % Magnification
Camera.f_num = 5.7; % F-number