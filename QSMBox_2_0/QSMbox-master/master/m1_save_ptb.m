function flagout = m1_save_ptb(ptb,flag)
% DO NOT EDIT
% Created by Julio Acosta-Cabronero

flagout = flag;

if strcmpi(flag,'preset')
    [~,ptb.settings_file]   = fileparts(which(mfilename));
    ptb.settings_file       = [ptb.settings_file '.m'];
    save('ptb.mat','ptb')

elseif strcmpi(flag,'run') || strcmpi(flag,'saveptb')
    save('ptb.mat','ptb')

elseif strcmpi(flag,'runextptb')
    flagout = 'run';

end
