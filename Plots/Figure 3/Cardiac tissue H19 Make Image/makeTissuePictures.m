clear

dp = readmm('dapi009.tif');
al = readmm('alexa009.tif');
cy = readmm('cy009.tif');
tm = readmm('tmr009.tif');

dp = max(dp.imagedata(:,:,5:10),[],3);
al = max(al.imagedata(:,:,:),[],3);
cy = max(cy.imagedata(:,:,:),[],3);
tm = max(tm.imagedata(:,:,:),[],3);

R = [107.9361   97.0246  492  610];

dp = imcrop(dp,R);
al = imcrop(al,R);
cy = imcrop(cy,R);
tm = imcrop(tm,R);

dp = scale(dp);
al = scale(al)*3.5;
cy = scale(cy)*1.5;
tm = scale(tm)*2.5;


cyan = [0 174 239]/256*1.5;  % *2 is to increase the brightness a bit
orange = [247 148 30]/256*1.5;

al2 = cat(3,al,al,al+dp*.75);
cy2 = cat(3,cy*cyan(1),cy*cyan(2),cy*cyan(3)+dp*.75);
tm2 = cat(3,tm*orange(1),tm*orange(2),tm*orange(3)+dp*.75);
dp2 = cat(3,0 *dp, 0* dp, dp*1.2);

imwrite(im2uint16(al2),'outGuide.tif');
imwrite(im2uint16(cy2),'outCy.tif');
imwrite(im2uint16(tm2),'outTm.tif');
imwrite(im2uint16(dp2),'outdp.tif');

R = [245.6542    9.5843  210 120];

imwrite(imcrop(im2uint16(al2),R),'outGuideZoom.tif');
imwrite(imcrop(im2uint16(cy2),R),'outCyZoom.tif');
imwrite(imcrop(im2uint16(tm2),R),'outTmZoom.tif');
imwrite(imcrop(im2uint16(dp2),R),'outdpZoom.tif');
