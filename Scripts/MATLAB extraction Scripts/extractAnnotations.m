outDir = '/Volumes/SNPFISH_H19/H19 SNP-FISH Paper/Raw Data/5-AZA/extractedCSV';

%%

expMetaData = ReadYaml('readMe.yaml');


dataExtractor = improc2.launchDataExtractor;
dataExtractor.extractFromProcessorData('alexa.RNACounts', @getNumSpots, 'alexa:Spots');

tools = improc2.launchImageObjectTools

objH = tools.objectHandle;
    getObjArea = @(objH) sum(sum(objH.getCroppedMask()));
    dataExtractor.extractFromObj('area', getObjArea);

cd(outDir)
dataExtractor.extractAllToCSVFile([expMetaData.name, '.csv'])
 