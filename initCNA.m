function [ returnValue ] = initCNA
%INITCNA Initializing CellNetAnalyzer
%   This function calls the init function of CellNetAnalyzer

directoriesToSearch=char( ...
    'W:/MetabolicModeling/CellNetAnalyzer/CellNetAnalyzer_for_MATLAB75_or_higher', ...
    'C:/Users/centlerf/Nextcloud2/Home/Projects/MetabolicModeling/FBATools/CellNetAnalyzer/CellNetAnalyzer_for_MATLAB75_or_higher', ...
    'U:/DBFZ/Florian_Kopplung/FBACoupling/CNA2012.1_for_MATLAB75_or_higher/CellNetAnalyzer_for_MATLAB75_or_higher' ...
   );

for i = 1:length(directoriesToSearch)
    if (exist(directoriesToSearch(i,:), 'dir')) 
        cwd = pwd();
        cd(directoriesToSearch(i,:));
        startcna(1);
        cd(cwd);
        returnValue = 0;
        return
    end
end

% failure to initialize CNA
returnValue = 1;

end

