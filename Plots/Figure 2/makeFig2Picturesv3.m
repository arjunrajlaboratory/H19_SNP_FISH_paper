clear

% First let's do the monoallelic one
cd('Transcript Raw Images_v3');
dp = readmm('dapi022.tif');
tm = readmm('tmr022.tif');
al = readmm('alexa022.tif');
cy = readmm('cy022.tif');
gfp = readmm('gfp022.tif');

% Just keep planes 9-16 for better clarity to the image.  Could also
% probably just keep a single plane if we want it to look cleaner.
dp1 = max(dp.imagedata(:,:,9:21),[],3);
tm1 = max(tm.imagedata(:,:,9:21),[],3);
al1 = max(al.imagedata(:,:,9:21),[],3);
cy1 = max(cy.imagedata(:,:,9:21),[],3);
gfp1 = max(gfp.imagedata(:,:,9:18),[],3);



R = [670.6693  345.5435  225  165];  % Crop region.
R = [55 280 225 175];

dp1 = imcrop(dp1,R); % Crop images.
tm1 = imcrop(tm1,R);
al1 = imcrop(al1,R);
cy1 = imcrop(cy1,R);
gfp1 = imcrop(gfp1,R);

% Scale all channels

dp1 = scale(dp1)/1;  % Tone down the DAPI.
al1 = scale(al1)*3;

gfp_forScale = readmm('gfp028.tif');
gfp_forScale = max(gfp_forScale.imagedata(:,:,9:18),[],3);
R_forScale = [650 545 225 175];
gfpforScale = imcrop(gfp_forScale,R);

cd ..

mngfp = min(gfp1(:));
maxgfp = max(gfp1(:));

gfp1 = scale(gfp1, 'intensityScale', [mngfp maxgfp]); %There is only about 50 pixel intensity spread in raw images

mn = min([tm1(:);cy1(:)]);
mx = max([tm1(:);cy1(:)]) - 200;  % Slightly oversaturate the high end so we can see individual spots a bit better.

tm1 = scale(tm1,'intensityScale',[mn mx]);  % Pick same scales for both TMR and CY for a better apples to apples comparison
cy1 = scale(cy1,'intensityScale',[mn mx]);

cyan = [0 174 239]/256*2;  % *2 is to increase the brightness a bit
orange = [247 148 30]/256*2;
green = [0 50 0]/256;

RGB = cat(3,cy1*orange(1),cy1*orange(2),cy1*orange(3)+dp1);
imwrite(im2uint16(RGB),'IGF2Bi_cyCASTv3-1.tif')

RGB = cat(3,al1,al1,al1+dp1);
imwrite(im2uint16(RGB),'IGF2Bi_AlexaGuidev3-1.tif')

RGB = cat(3, tm1*cyan(1),tm1*cyan(2),tm1*cyan(3)+dp1);
imwrite(im2uint16(RGB),'IGF2Bi_tmrB6v3-1.tif')

RGB = cat(3, gfp1 + green(1), gfp1 + green(2), gfp1 + green(3) + dp1);
imwrite(im2uint16(RGB),'IGF2Bi_gfpIGF2v3-1.tif')
%%

clear

% Second Version
cd('Transcript Raw Images_v3');
dp = readmm('dapi028.tif');
tm = readmm('tmr028.tif');
al = readmm('alexa028.tif');
cy = readmm('cy028.tif');
gfp = readmm('gfp028.tif');

% Just keep planes 9-16 for better clarity to the image.  Could also
% probably just keep a single plane if we want it to look cleaner.
dp1 = max(dp.imagedata(:,:,9:18),[],3);
tm1 = max(tm.imagedata(:,:,9:18),[],3);
al1 = max(al.imagedata(:,:,9:18),[],3);
cy1 = max(cy.imagedata(:,:,9:18),[],3);
gfp1 = max(gfp.imagedata(:,:,9:18),[],3);


cd ..

R = [670.6693  345.5435  225  165];  % Crop region.
R = [650 545 225 175];

dp1 = imcrop(dp1,R); % Crop images.
tm1 = imcrop(tm1,R);
al1 = imcrop(al1,R);
cy1 = imcrop(cy1,R);
gfp1 = imcrop(gfp1,R);

% Scale all channels

dp1 = scale(dp1)/1;  % Tone down the DAPI.
al1 = scale(al1)*4.5;
gfp1 = scale(gfp1) * 3;

mn = min([tm1(:);cy1(:)]);
mx = max([tm1(:);cy1(:)]) - 400;  % Slightly oversaturate the high end so we can see individual spots a bit better.

tm1 = scale(tm1,'intensityScale',[mn mx]);  % Pick same scales for both TMR and CY for a better apples to apples comparison
cy1 = scale(cy1,'intensityScale',[mn mx]);

cyan = [0 174 239]/256*2;  % *2 is to increase the brightness a bit
orange = [247 148 30]/256*2;
green = [0 50 0]/256;

RGB = cat(3,cy1*orange(1),cy1*orange(2),cy1*orange(3)+dp1);
imwrite(im2uint16(RGB),'IGF2mono_cyCASTv3-1.tif')

RGB = cat(3,al1,al1,al1+dp1);
imwrite(im2uint16(RGB),'IGF2mono_AlexaGuidev3-1.tif')

RGB = cat(3, tm1*cyan(1),tm1*cyan(2),tm1*cyan(3)+dp1);
imwrite(im2uint16(RGB),'IGF2monof_tmrB6v3-1.tif')

RGB = cat(3, gfp1 + green(1), gfp1 + green(2), gfp1 + green(3) + dp1);
imwrite(im2uint16(RGB),'IGF2monof_gfpIGF2v3-1.tif')

