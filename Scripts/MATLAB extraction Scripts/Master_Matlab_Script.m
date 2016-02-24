% Matlab Master Script for H19 SNP-FISH Paper. Run to reproduce, from raw
% data, all CSV files that are subsequently imported into R for plot
% generation. 
%
%




%Make sure that this is the address of the H19 SNP-FISH Folder on your
%local computer
currentDir = '/Volumes/Paul_Travel/H19 Paper Submssion_ALL_RAW_Images/H19 SNP-FISH Paper';



%% Extract 5-AZA Count Data

currentDirAZA = [currentDir, '/Raw Data/5-AZA'];
outDir = [currentDirAZA, '/extractedCSV'];
cd(currentDirAZA)


rawDataContents = struct2dataset(dir([currentDirAZA]));
rawDataDirNames = cellstr(rawDataContents.name);

rawDataFolderIndex =  cell2mat(cellfun(@(x) ~isempty(regexp(x, '^Rep')), ...
            rawDataDirNames,'UniformOutput', false));

rawDataFolderNames = rawDataDirNames(rawDataFolderIndex);

for i = 1:numel(rawDataFolderNames)
    tempdir = [currentDirAZA, '/', rawDataFolderNames{i}];
    cd(tempdir)
    rawDatasubContents = struct2dataset(dir(pwd));
    rawDatasubDirNames = cellstr(rawDatasubContents.name);
    
    rawDatasubFolderIndex =  cell2mat(cellfun(@(x) isempty(regexp(x, '^\.')), ...
        rawDatasubDirNames,'UniformOutput', false));
    
    rawDatasubFolderNames = rawDatasubDirNames(rawDatasubFolderIndex);
    
    for j = 1:numel(rawDatasubFolderNames)
        cd([tempdir, '/', rawDatasubFolderNames{j}])
        
        expMetaData = ReadYaml('readMe.yaml');
        
        
        dataExtractor = improc2.launchDataExtractor;
        dataExtractor.extractFromProcessorData('alexa.RNACounts', @getNumSpots, 'alexa:Spots');
        
        tools = improc2.launchImageObjectTools;
        
        objH = tools.objectHandle;
        getObjArea = @(objH) sum(sum(objH.getCroppedMask()));
        dataExtractor.extractFromObj('area', getObjArea);
        
        cd(outDir)
        dataExtractor.extractAllToCSVFile([expMetaData.name, '.csv'])
    end
end


%% Extact Bisulfite Clone Data
 

IGF2CloneDir = [currentDir, '/Raw Data/Clone Bisulfite/07-27-2015_IGF2+_1-3_4x'];
BiCloneDir = [currentDir, '/Raw Data/Clone Bisulfite/06-04-2015_Bi_5-3_1x_Bi#1_Fig4'];

outCloneDir = [currentDir, '/Raw Data/Clone Bisulfite/ExtractedCSVs'];

cd(IGF2CloneDir)
extractSpotCSVs_gfp_cellArrayObj(outCloneDir)
cd(currentDir)

cd(BiCloneDir)
extractSpotCSVs_gfp_cellArrayObj(outCloneDir)
cd(currentDir)


%% Extract raw Data Folders names in Directory for all other Data
cd([currentDir, '/Raw Data'])
outDir = [currentDir, '/Raw Data/ExtractedCSVs'];


rawDataContents = struct2dataset(dir([currentDir, '/Raw Data']));
rawDataDirNames = cellstr(rawDataContents.name);

%All data folders start with "Data_"
rawDataFolderIndex =  cell2mat(cellfun(@(x) ~isempty(regexp(x, '^Data\_')), ...
            rawDataDirNames,'UniformOutput', false));

rawDataFolderNames = rawDataDirNames(rawDataFolderIndex);



IGF2DataFolderIndex = cell2mat(cellfun(@(x) ~isempty(regexp(x, 'IGF2')), ...
            rawDataFolderNames,'UniformOutput', false));
        
IGF2DataFolderNames = rawDataFolderNames(IGF2DataFolderIndex);
NONIGF2DataFolderNames = rawDataFolderNames(~IGF2DataFolderIndex);
%% Extract IGF2 DATA
for i = 1:numel(IGF2DataFolderNames)
    tempdir = [currentDir, '/Raw Data/', IGF2DataFolderNames{i}];
    cd(tempdir)
    extractSpotCSVs_gfp_cellArrayObj(outDir) 
    cd(currentDir)
end

%% Extract NON IGF2 DATA

for i = 1:numel(NONIGF2DataFolderNames)


tempdir = [currentDir, '/Raw Data/', NONIGF2DataFolderNames{i}];
 
cd(tempdir)
extractSpotCSVs_cellArrayObj(outDir)
cd(currentDir)
end


%% Extract Single Oligo TEST

DataDir = '/Users/paulginart/Dropbox/H19 SNP-FISH Paper/Raw Data/Single Oligo Test MEF 9CG';
outDir = '/Users/paulginart/Dropbox/H19 SNP-FISH Paper/Raw Data/ExtractedCSVs';

rawDataContents = struct2dataset(dir(DataDir));
rawDataDirNames = cellstr(rawDataContents.name);

rawDataFolderIndex =  cell2mat(cellfun(@(x) ~isempty(regexp(x, '^9CG')), ...
            rawDataDirNames,'UniformOutput', false));

rawDataFolderNames = rawDataDirNames(rawDataFolderIndex);

for i = 1:numel(rawDataFolderNames)
    tempdir = [DataDir, '/', rawDataFolderNames{i}];
    cd(tempdir)
    extractSpotCSVs_gfp_cellArrayObj(outDir)
    cd(currentDir)
end


%% Extract Figure2 WT

DataDir = '/Volumes/SNPFISH_H19/H19 SNP-FISH Paper/Raw Data/Data_CxB_MEF_WT_IGF2_2';
outDir = '/Volumes/SNPFISH_H19/H19 SNP-FISH Paper/Raw Data/ExtractedCSVs';

cd(DataDir)
extractSpotCSVs_gfp_cellArrayObj(outDir)
cd(currentDir)

%% Extract Figure 3 WT



DataDir = '/Volumes/SNPFISH_H19/H19 SNP-FISH Paper/Raw Data/Data_CxB_9CG_MEFinCRL';
outDir = '/Volumes/SNPFISH_H19/H19 SNP-FISH Paper/Raw Data/ExtractedCSVs';

cd(DataDir)
extractSpotCSVs_cellArrayObj(outDir)
%cd(currentDir)



