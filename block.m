
function [block_max, block_min]=block(I)

maxfun = @(block_struct) max(block_struct.data);

block_size = [50 50];
block_max = blockproc(I,block_size,maxfun);
disp(block_max);
figure
imshow(block_max);
title('Block Processing - Simplest Syntax');

minfun = @(block_struct) min(block_struct.data);

block_size = [50 50];
block_min = blockproc(I,block_size,minfun);
disp(block_min);
figure
imshow(block_min);
title('Block Processing - Simplest Syntax');