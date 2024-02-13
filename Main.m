clear
close all
currPath = fileparts(mfilename('fullpath'));% get current path
p = cd(currPath);  %open the m file folder
addpath(genpath( './Main_fun'));
my_pro;