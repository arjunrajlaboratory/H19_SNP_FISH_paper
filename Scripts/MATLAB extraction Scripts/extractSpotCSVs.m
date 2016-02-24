function [] = extractSpotCSVs(outDir)
%Extracts Spot CSVs after SNP_FISH_Processing
%   outDir- data output directory

expMetaData = ReadYaml('readMe.yaml'); %Read in YAML file for experiment details


dataFiles = improc2.utils.ImageObjectDataFiles();
onDisk = improc2.utils.FileBasedImageObjectArrayCollection(dataFiles);

tools = improc2.launchImageObjectTools(onDisk);


snpMap.channels = {'alexa', 'tmr', 'cy'};
snpMap.names = {expMetaData.channels.alexa.probe, ...
                 expMetaData.channels.tmr.probe, ...
                 expMetaData.channels.cy.probe,};
             
expName = expMetaData.name;
             
SNP_Processor = expMetaData.SNP_Processor;
PS_Processor = expMetaData.PS_Processor;

% Create Data Tables for this Experiment

cellCounts = [];
tools.iterator.goToFirstObject
cellID = 1;

alexaAll = []; 
cyAll = [];
tmrAll = [];

alexa_PS = []; 
cy_PS = [];
tmr_PS = [];



while(tools.iterator.continueIteration)
    results = tools.objectHandle.getData(SNP_Processor);
 
    
    alexa = struct2table(results.data.alexa);
    alexa.cellID = ones(height(alexa),1) * cellID;
    cy = struct2table(results.data.cy);
    cy.cellID = ones(height(cy),1) * cellID;
    tmr = struct2table(results.data.tmr);
    tmr.cellID = ones(height(tmr),1) * cellID;
    
    if isempty(alexaAll)
    alexaAll = alexa;
    else
    alexaAll = vertcat(alexaAll, alexa);
    end
    
    if isempty(cyAll)
    cyAll = cy;
    else
    cyAll = vertcat(cyAll, cy);
    end
    
    if isempty(tmrAll)
        tmrAll = tmr; 
    else
    tmrAll = vertcat(tmrAll, tmr);
    end
    
    
    results_PS = tools.objectHandle.getData(PS_Processor);
    
    alexa2 = struct2table(results_PS.data.alexa);
    alexa2.cellID = ones(height(alexa2),1) * cellID;
    cy2 = struct2table(results_PS.data.cy);
    cy2.cellID = ones(height(cy2),1) * cellID;
    tmr2 = struct2table(results_PS.data.tmr);
    tmr2.cellID = ones(height(tmr2),1) * cellID;
    
    if isempty(alexa_PS)
    alexa_PS = alexa2;
    else
    alexa_PS = vertcat(alexa_PS, alexa2);
    end
    
    if isempty(cy_PS)
    cy_PS = cy2;
    else
    cy_PS = vertcat(cy_PS, cy2);
    end
    
    if isempty(tmr_PS)
        tmr_PS = tmr2; 
    else
    tmr_PS = vertcat(tmr_PS, tmr2);
    end
    
    tools.iterator.goToNextObject;
    cellID = cellID + 1;
    
end



% Update Table Names
currentNames = alexaAll.Properties.VariableNames;
updateNames = currentNames;

updateNames{6} = [snpMap.channels{2}, '_ID'];
updateNames{7} = [snpMap.channels{2}, '_positions'];
updateNames{8} = [snpMap.channels{2}, '_amplitude'];
updateNames{9} = [snpMap.channels{2}, '_sigma'];
updateNames{10} = [snpMap.channels{3}, '_ID'];
updateNames{11} = [snpMap.channels{3}, '_poistions'];
updateNames{12} = [snpMap.channels{3}, '_amplitude'];
updateNames{13} = [snpMap.channels{3}, '_sigma'];

alexaAll.Properties.VariableNames = updateNames;
alexa_PS.Properties.VariableNames = updateNames;


cd(outDir)
writetable(alexaAll,[expName '_MasterGuideSpotsTable.csv'])
writetable(alexa_PS,[expName '_MasterGuideSpotsTable_Shifted.csv'])


end

