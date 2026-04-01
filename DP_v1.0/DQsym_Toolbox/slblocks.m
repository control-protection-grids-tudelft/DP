function blkStruct = slblocks
% Add DQN_Lib to the Simulink Library Browser

Browser.Library = 'DQsym_Lib';      % name of the SLX file without extension
Browser.Name    = 'DQsym Library';  % name shown in Library Browser
Browser.IsFlat  = 0;              % 0 if library has subsystems/categories

blkStruct.Browser = Browser;
end