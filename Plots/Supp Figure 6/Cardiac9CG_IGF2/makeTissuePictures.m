clear

dp = readmm('dapi001.tif');
al = readmm('alexa001.tif');
cy = readmm('cy001.tif');
tm = readmm('tmr001.tif');
gf = readmm('gfp001.tif');

dp = max(dp.imagedata(:,:,13:29),[],3);
al = max(al.imagedata(:,:,13:29),[],3);
cy = max(cy.imagedata(:,:,13:29),[],3);
tm = max(tm.imagedata(:,:,13:29),[],3);
gf = max(gf.imagedata(:,:,13:29),[],3);

R = [195  50  175  580];

dp = imcrop(dp,R);
al = imcrop(al,R);
cy = imcrop(cy,R);
tm = imcrop(tm,R);
gf = imcrop(gf,R);

dp = scale(dp);
al = scale(al)*2.5;
%cy = scale(cy)*1.1;
%tm = scale(tm)*1.1;
gf = scale(gf)*1;

mn = min([tm(:);cy(:)]);
mx = max([tm(:);cy(:)]) - 200;  % Slightly oversaturate the high end so we can see individual spots a bit better.

tm = scale(tm,'intensityScale',[mn mx-300]);  % Pick same scales for both TMR and CY for a better apples to apples comparison
cy = scale(cy,'intensityScale',[mn mx]);

mn = min([tm(:);cy(:)]);
mx 


cyan = [0 174 239]/256*2.1;  % *2 is to increase the brightness a bit
orange = [247 148 30]/256*1.5;
green = [0 80 0]/256;

al2 = cat(3,al,al,al+dp*.75);
cy2 = cat(3,cy*orange(1),cy*orange(2),cy*orange(3)+dp*.75);
tm2 = cat(3,tm*cyan(1),tm*cyan(2),tm*cyan(3)+dp*.75);
gf2 = cat(3, gf + green(1), gf + green(2), gf + green(3) + dp);

imwrite(im2uint16(al2),'outGuide.tif');
imwrite(im2uint16(cy2),'outCy.tif');
imwrite(im2uint16(tm2),'outTm.tif');
imwrite(im2uint16(gf2),'outGfp.tif');

