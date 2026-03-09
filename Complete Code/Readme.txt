Three-Dimensional Dynamic Speckle Simulation Toolbox

File Structure
Core Files：
1.ModelSetting.m - Model parameter configuration file
2.dEMC.m - Main dynamic Electric field Monte Carlo simulation 
3.demo.m - Parallel running and OU process

Mie Scattering Calculation Module:
mie.m - Mie scattering efficiency calculation
mie_abcd.m - Mie coefficients calculation
Mie_pt.m - Angular function calculation
Mie_S12.m - Mie scattering amplitude calculation
PreparMieScattering.m - Mie scattering parameter preprocessing
tissueTypeMapping.m - Tissue type mapping function

System Requirements：
1.MATLAB R2019a or later
2.Parallel Computing Toolbox (for parallel computation)
3.Multi-core CPU recommended for better performance

Installation
Download all MATLAB files into a single directory
Add the directory to MATLAB's search path

Basic Usage Workflow
Step 1: Configure Simulation Parameters
Edit ModelSetting.m to set your simulation parameters.

Step 2: Run demo.m in parallel


Primary Outputs:
MieC: Contains Mie scattering coefficients
  ├── n_sctCOMPLEX  % Scatterer refractive indices
  ├── dia           % Scatterer diameters
  ├── lsa           % Scattering mean free paths
  ├── lt            % Transport mean free paths
  └── gFactor       % Anisotropy factors

iImage: Contains 3D image data (per simulation run)
  ├── IMRS_layer_pol  % Single scattering reflection (polarized)
  ├── IMRM_layer_pol  % Multiple scattering reflection (polarized)
  ├── G1_perm_RS_layer % Reflection g1 values
  ├── Photonpath_perm_RS_layer % Reflection photon counts
  └── ... (16 total variables)

Image: Combined results from multiple runs (output of demo.m)
  ├── All combined image fields from iImage
  ├── G1RS_layer    % Normalized reflection g1 (single scattering)
  ├── G1RM_layer    % Normalized reflection g1 (multiple scattering)
  └── ... additional normalized fields

RIs    %Reflection dynamic speckle for single scattering
RIm    %Reflection dynamic speckle for multiple scattering
TIs    %Transmission dynamic speckle for single scattering
TIm    %Transmission dynamic speckle for multiple scattering