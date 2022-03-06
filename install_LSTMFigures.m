function install_LSTMFigures()
[filepath,name,ext] = fileparts(mfilename('fullpath'));
addpath(genpath(filepath));
if(isempty(which('syspec_ssim_v15')))
    error('All the necessary functions have not been successfully added to the MATLAB path.');
else
    disp('LSTMFigures package has been successfully installed.');
end
end