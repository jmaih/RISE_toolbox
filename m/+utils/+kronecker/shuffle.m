%--- help for matlab.io.datastore.ImageDatastore/shuffle ---
%
% shuffle Shuffle the files of a copy of the ImageDatastore.
%    DSOUT = SHUFFLE(DS) creates a deep copy of the input ImageDatastore and
%    shuffles the files using randperm, resulting in the datastore DSOUT.
% 
%    Example:
%    --------
%       folders = fullfile(matlabroot,'toolbox','matlab',{'demos','imagesci'});
%       exts = {'.jpg','.png','.tif'};
%       imds = imageDatastore(folders,'LabelSource','foldernames','FileExtensions',exts)
%       shuffledDs = shuffle(imds)
% 
%    See also imageDatastore, splitEachLabel, countEachLabel, hasdata,
%    readimage, readall, preview, reset.
%
%    Documentation for matlab.io.datastore.ImageDatastore/shuffle
%       doc matlab.io.datastore.ImageDatastore/shuffle
%
%    Other uses of shuffle
%
%       matlab.io.datastore.Shuffleable/shuffle
%