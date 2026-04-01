function install_DQN_Lib()
% Install DQN_Lib into the Simulink Library Browser

rootFolder = fileparts(mfilename('fullpath'));

% Add the toolbox folder to path
addpath(rootFolder);

% Save path for future MATLAB sessions
savepath;

% Refresh Simulink customizations and Library Browser
rehash toolboxcache;
sl_refresh_customizations;

% Refresh Library Browser window if available
try
    lb = LibraryBrowser.LibraryBrowser2;
    refresh(lb);
catch
end

% Open the library once to verify it is accessible
open_system('DQsym_Lib');

disp('DQsym_Lib installed.');
disp('Open the Simulink Library Browser and look for "DQN Library".');
end