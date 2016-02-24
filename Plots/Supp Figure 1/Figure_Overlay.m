%% Create Sample View Region for Figure 1 showing colocalization

clear;
improc2.tests.cleanupForTests;

expMetaData = ReadYaml('readMe.yaml'); %Read in YAML file for experiment details


dataFiles = improc2.utils.ImageObjectDataFiles();
onDisk = improc2.utils.FileBasedImageObjectArrayCollection(dataFiles);
collection = improc2.utils.loadCollectionIntoMemory(onDisk);

tools = improc2.launchImageObjectTools(onDisk);



snpMap.channels = {'alexa', 'tmr', 'cy'};
snpMap.names = {expMetaData.channels.alexa.probe, ...
                expMetaData.channels.tmr.probe, ...
                expMetaData.channels.cy.probe,};

            

%%


results = tools.objectHandle.getData('snpColoc');
guide = results.data.(snpMap.channels{1});

img_tmr = readmm('tmr003.tif');
img_Alexa = readmm('alexa003.tif');
img_cy = readmm('cy003.tif');
img_dapi = readmm('dapi003.tif');



R = [120 0 470 750];

imgmax_Alexa = max(img_Alexa.imagedata(:,:,7:12),[],3);  %MAX MERGE
imgmax_cy = max(img_cy.imagedata(:,:,7:12),[],3);  %MAX MERGE
imgmax_tmr = max(img_tmr.imagedata(:,:,7:12),[],3);  %MAX MERGE
imgmax_dapi = max(img_dapi.imagedata(:,:,7:12),[],3);  %MAX MERGE


cy1 = imcrop(imgmax_cy,R);
al1 = imcrop(imgmax_Alexa,R);
tm1 = imcrop(imgmax_tmr,R);
dp1 = imcrop(imgmax_dapi,R);

al1 = scale(al1)*1.9;
dp1 = scale(dp1);

mn = min([tm1(:);cy1(:)]);
mx = max([tm1(:);cy1(:)]);  

cyan = [0 174 239]/256*2;  % *2 is to increase the brightness a bit
orange = [247 148 30]/256*2;

tm1 = scale(tm1, [mn mx]);%,'intensityScale',[mn mx]);  % Pick same scales for both TMR and CY for a better apples to apples comparison
cy1 = scale(cy1, [mn mx]);

RGB = cat(3,cy1 * orange(1),cy1 * orange(2),cy1 * orange(3));
imwrite(im2uint16(RGB),'Fig1_Cy_C7.tif')

RGB = cat(3,al1,al1,al1 + dp1);
imwrite(im2uint16(RGB),'Fig1_AlexaGuide.tif')

RGB = cat(3,tm1 * cyan(1),tm1 * cyan(2),tm1 * cyan(3));
imwrite(im2uint16(RGB),'Fig1_tmr_B6.tif')

%%
fighandle = figure;

img = tools.objectHandle.getData([snpMap.channels{1}, ':Spots']).getImage;
imgmax = max(img,[],3);  %MAX MERGE
imgmax = imadjust(imgmax, stretchlim(imgmax, [0 0.995]));


figure(fighandle)
imshow(imgmax,[]);  %Plot image
hold on
plot(guide.position(:,1), guide.position(:,2), 'wo'); 
plot(guide.position(guide.labels == 'C7',1), guide.position(guide.labels == 'C7',2),'co','markersize',6); %Plot the co_localized spots
plot(guide.position(guide.labels ==  'B6',1), guide.position(guide.labels == 'B6',2),'Color', [1, 0.65, 0.2], ...
    'LineStyle', 'none', 'Marker','o','markersize',6); %Plot the co_localized spots
plot(guide.position(guide.labels == '3-color',1), guide.position(guide.labels == '3-color',2),'mo','markersize',6);
legend({snpMap.names{1}, snpMap.names{2},snpMap.names{3},'3-color'})
hold off
%%
print('-dtiffn', 'Fig1_Overlay.tif');

imgOverlay  = imread('Fig1_Overlay.tif');
imwrite(imgOverlay,'Fig1_Overlay.tif')

