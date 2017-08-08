function moments = BBC(image)


imgGray = double(rgb2gray(image))/255;
imgGray = imresize(imgGray, [256 256]);
imbin   = im2bw(image);
coeff_0 = dwt2(imbin', 'coif1');
coeff_1 = dwt2(imgGray', 'coif1');
coeff_2 = dwt2(coeff_1, 'coif1');
coeff_3 = dwt2(coeff_2, 'coif1');
coeff_4 = dwt2(coeff_3, 'coif1');

% construct the feaute vector
meanCoeff = mean(coeff_4);
stdCoeff = std(coeff_4);

moments = [meanCoeff stdCoeff];

end