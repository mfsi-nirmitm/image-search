function similarity(numOfReturnedImages, queryImageFeatureVector, dataset)

query_image_name = queryImageFeatureVector(:, end);
dataset_image_names = dataset(:, end);

queryImageFeatureVector(:, end) = [];
dataset(:, end) = [];

% compute manhattan distance
manhattan = zeros(size(dataset, 1), 1);
for k = 1:size(dataset, 1)

    manhattan(k) = sum( abs(dataset(k, :) - queryImageFeatureVector) ./ ( 1 + dataset(k, :) + queryImageFeatureVector ) );
end


manhattan = [manhattan dataset_image_names];


[sortedDist indx] = sortrows(manhattan);
sortedImgs = sortedDist(:, 2);


arrayfun(@cla, findall(0, 'type', 'axes'));

% display similarity image
str_name = int2str(query_image_name);
queryImage = imread( strcat('images\', str_name, '.jpg') );
subplot(3, 7, 1);
imshow(queryImage, []);
title('Query Image', 'Color', [1 0 0]);

% dispaly images returned by similarity
for m = 1:numOfReturnedImages
    img_name = sortedImgs(m);
    img_name = int2str(img_name);
    str_name = strcat('images\', img_name, '.jpg');
    returnedImage = imread(str_name);
    subplot(3, 7, m+1);
    imshow(returnedImage, []);
end

end