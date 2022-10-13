% script2_find finds the those parameter values that were searched from the
% script1_dosage file and sorts those values that show certain levels of
% sog and sna borders and also meet the hard constraint that the
% concentration of dl0 = 1 at dorsal border be greater that 0.5 at the
% final time point.

clear
clc
close all

%++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
% CONTROLS
Results_folder = '../Results/Runs4/3/';         % results are in folder 1 and 3
% Results_folder = 'Results/';
error           = 1.5;
hard_constraint = 0.5;


Mats_folder  = strcat(Results_folder,'Mats/');
Figs_folder  = strcat(Results_folder,'Figs/');

%++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++


%------------------------------------
% Extracting folders and Initializing
matfiles     = dir(fullfile(Mats_folder));
nfolders     = length(matfiles);   
Found_Params = [];
AllParamsTrials    = [];
Found_Errors = [];
Found_Concs  = [];
count_found    = 1;
count_unfound  = 1;
total        = 0;


%-------------------------
% looping over all folders
% NOTE - Starts from 3 to remove '.' and '..'
for i = 3:nfolders 
   
    flag = 1;
    %-----------------------
    % Read files from folder
    filesdir = strcat(Mats_folder, matfiles(i).name,'/');
    try
        load(strcat(filesdir,'AllParams.mat'),'AllParams');
        load(strcat(filesdir,'AllErrors.mat'),'AllErrors');
        load(strcat(filesdir,'AllConcs.mat'),'AllConcs');
        AllParamsTrials = [AllParamsTrials; AllParams]; %#ok<AGROW>
        total = total + length(AllParams);
    catch
        fprintf('Params not found in %i\n', i);
        flag = 0;
    end
    
    
      
    %---------------------------
    % Check if error values work
    if flag == 1
        n     = length(AllErrors);
        for j = 1:n
            if (AllErrors(j,1:3) < error)  
                if (AllErrors (j,4) > hard_constraint)
                    if AllConcs(j,:)>= 0
                        Found_Params(count_found,:) = AllParams(j,:);
                        Found_Errors(count_found,:) = AllErrors(j,:);
                        Found_Concs(count_found,:)  = AllConcs(j,:);
                        count_found = count_found + 1;
                    end
                end
            else
                unFound_Params(count_unfound,:) = AllParams(j,:);
                unFound_Errors(count_unfound,:) = AllErrors(j,:);
                unFound_Concs(count_unfound,:)  = AllConcs(j,:);
                count_unfound = count_unfound + 1;
            end         
        end
    end
end

% ------------------------------------
% Save values that satisfy conditions!
save(strcat(Results_folder,'Found_Params.mat'),'Found_Params');
save(strcat(Results_folder,'Found_Errors.mat'),'Found_Errors');
save(strcat(Results_folder,'Found_Concs.mat'),'Found_Concs');

save(strcat(Results_folder,'unFound_Params.mat'),'unFound_Params');
save(strcat(Results_folder,'unFound_Errors.mat'),'unFound_Errors');
save(strcat(Results_folder,'unFound_Concs.mat'),'unFound_Concs');

save(strcat(Results_folder,'AllParamsTrials.mat'),'AllParamsTrials');

