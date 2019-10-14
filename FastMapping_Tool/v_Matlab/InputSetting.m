%% Inputs file for FastMapping_Tool

% Example iclude a visual stimuli task (single-voxel) 





% Data Inputs

func_name = 'Visual_Voxel';    % Visual_Voxel
func_file = fullfile('DataExample/',['TCN_',func_name,'.mat']);

load(func_file);

% Data Info


param.Tpoints=size(TCN,1);
param.NbrVoxels = size(TCN,2);

if (param.Tpoints > param.NbrVoxels)
  fprintf(['Number of time points is larger than number of voxels ! \n',...
  'Is this ok ? Otherwise, maybe you want to take the transpose of the input matrix. \n'])
end


 % %%%%%%%%%%%% SET BY USER %%%%%%%%%%%%%%
param.TR = 2.4;        



% parameter setting

param.niter = 400;               % Fista iterations
param.condition = 'blocks';     % OR condition = 'spikes';
param.lam0 = 1;                 % weight for temporal regularization, the internal iteration 
                                % uses wavelet noise estimation
                                
                                
                                 



