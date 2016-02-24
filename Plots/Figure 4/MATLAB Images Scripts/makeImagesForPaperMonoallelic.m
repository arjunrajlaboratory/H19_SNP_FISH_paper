%% Read the data files
cd 07-27-2015_IGF2+_1-3_4x_ColorSwap
cy    = readmm('cy003.tif');
tmr   = readmm('tmr003.tif');
dapi  = readmm('dapi003.tif');
alexa = readmm('alexa003.tif');
gfp   = readmm('gfp003.tif');

%% Just keep the image data and drop all the metadata in the image file
cy    = cy.imagedata;
tmr   = tmr.imagedata;
dapi  = dapi.imagedata;
alexa = alexa.imagedata;
gfp   = gfp.imagedata;

%% Take the max-merges
mxcy = max(cy(:,:,5:10),[],3);
% Note: if you want just a few planes (e.g., 10-20), use the following:
%mxcy = max(cy(:,:,10:20),[],3);
mxtmr = max(tmr(:,:,5:10),[],3);
mxalexa = max(alexa(:,:,5:10),[],3);
mxgfp = max(gfp,[],3);
% For DAPI, let's just take a particular slice
mxdapi = dapi(:,:,8);

%% Use the following lines to get the rectangle of interest
% imshow(mxcy,[])
% R = getrect;
% R % This shows what the rectangle is

%% Now let's crop the images
% Once you've selected R, you can just use the following.
% I clean up the last two numbers to something relatively nice, because
% those are the width and height of the cropped region
R = [354.9615  607.9779  286.4336  281.6597];
R = [257.0967  469.5350  350 350];

cropCy = imcrop(mxcy,R);
cropTmr = imcrop(mxtmr,R);
cropAlexa = imcrop(mxalexa,R);
cropGfp = imcrop(mxgfp,R);
cropDapi = imcrop(mxdapi,R);

% Use statements like the following to show the Cy image
imshow(cropCy,[]);
% The [] part tells imshow to autoscale things.



%% Now time to scale the images
cropCyScaled = scale(cropCy);
cropTmrScaled = scale(cropTmr);
cropAlexaScaled = scale(cropAlexa);
cropGfpScaled = scale(cropGfp);
cropDapiScaled = scale(cropDapi);

% Now we don't need the [] in imshow:
imshow(cropCyScaled);
% If I wanted to pep things up a bit:
imshow(cropCyScaled*1.5);


%% Let's convert the images to doubles.
% This will make them easier to adjust.

cropCyScaled = im2double(cropCyScaled);
cropTmrScaled = im2double(cropTmrScaled)*2; % Factor of two to take care of max intensity differences.
cropAlexaScaled = im2double(cropAlexaScaled)*2;
cropGfpScaled = im2double(cropGfpScaled);
cropDapiScaled = im2double(cropDapiScaled);

%% Now let's make the color merge

% Each component is Red, Green, then Blue.
% Add to each component for white.
RGB = cat(3,cropCyScaled, cropCyScaled, cropCyScaled + cropDapiScaled);
imshow(RGB);

% Let's say we want to turn down the DAPI.  Just multiply it:
%RGB = cat(3,cropCyScaled, cropCyScaled, cropCyScaled + cropDapiScaled*0.5);
%imshow(RGB);

%% Now that we have a pic we like, let's write to disk:


cyan = [0 174 239]/256;  % *2 is to increase the brightness a bit
orange = [247 148 30]/256;

cd ..
% First, convert to 16 bit
RGB = cat(3,cropCyScaled*cyan(1), cropCyScaled*cyan(2), cropCyScaled*cyan(3) + cropDapiScaled);
imshow(RGB);
RGB = im2uint16(RGB);
imwrite(RGB,'monoalleliccyImage.tif');

RGB = cat(3,cropTmrScaled*orange(1), cropTmrScaled*orange(2), cropTmrScaled*orange(3) + cropDapiScaled);
imshow(RGB)
RGB = im2uint16(RGB);
imwrite(RGB,'monoallelictmrImage.tif');

RGB = cat(3,cropAlexaScaled, cropAlexaScaled, cropAlexaScaled + cropDapiScaled);
imshow(RGB)
RGB = im2uint16(RGB);
imwrite(RGB,'monoallelicalexaImage.tif');

