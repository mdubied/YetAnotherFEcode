clc

if ispc
    sslash = '\';
elseif isunix
    sslash = '/';
end

addpath(genpath(strcat(pwd,sslash,'src')));
addpath(genpath(strcat(pwd,sslash,'external')));
addpath(genpath(strcat(pwd,sslash,'examples',sslash,'Meshes')));

fprintf('\n------------------------------------------')
fprintf('\n|   *** Welcome to YetAnotherFEcode ***  |')
fprintf('\n------------------------------------------\n\n')