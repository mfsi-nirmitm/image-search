function [Y,Combmax,Combmin] = ODBTC(I)
  

red = I(:,:,1); % Red channel
green = I(:,:,2); % Green channel
blue = I(:,:,3); % Blue channel
a = zeros(size(I, 1), size(I, 2));
just_red = cat(3, red, a, a);
just_green = cat(3, a, green, a);
just_blue = cat(3, a, a, blue);
back_to_original_img = cat(3, red, green, blue);
figure, imshow(I), title('Original image')
figure, imshow(just_red), title('Red channel')
figure, imshow(just_green), title('Green channel')
figure, imshow(just_blue), title('Blue channel')
figure, imshow(back_to_original_img), title('Back to original image')

[Rblock_max, Rblock_min]=block(just_red);
[Bblock_max, Bblock_min]=block(just_green);
[Gblock_max, Gblock_min]=block(just_blue);

Combmax = cat(3, Rblock_max, Gblock_max, Bblock_max);
Combmin = cat(3, Rblock_min, Gblock_min, Bblock_min);

    % make sure I is double
    Iband = im2double(I);
 Ibar = mean(Iband, 3); %interband average image
 figure, imshow(Ibar), title('interband average image')
 
 if size(I,3)==3
    I=rgb2gray(I);
end
%size of image
[M,N]=size(I);
%convert to double
I=double(I);
Y=zeros(M,N);
%% Encoding
blksize=2;    %Block Size
mu=colfilt(I,[blksize,blksize],'distinct',@(x) ones(blksize^2,1)*mean(x));
sigma=colfilt(I,[blksize,blksize],'distinct',@(x) ones(blksize^2,1)*std(x));
q=I>mu;
q=colfilt(q,[blksize,blksize],'distinct',@(x) ones(blksize^2,1)*sum(x));
m=blksize^2;                          %length*width of block
a=mu-sigma.*(sqrt(q./m-q));           %low mean
b=mu+sigma.*(sqrt(m-q./q));           %high mean
H=I>=mu;                              %elements of Bitmap
Y(H)=a(H);
Y(~H)=b(~H);
Y=uint8(Y);                           %output BTC image
figure,imshow(Y)

 
 
 
 