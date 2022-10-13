function [IM,lsminf1,lsminf2] = dosage_imread(filename)
%
%
%
%
% This function reads in the images of the dosage experiment
%
% This function has been moved into the "Functions" folder.


lsm =  lsmRead2(filename);
lsminf1 = lsm.lsm;
lsminf2 = lsminfo(filename);
IM = lsm.data;
