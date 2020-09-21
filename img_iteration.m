files = dir('*.jpg');
var=0;
x=0;
a=false;
% s'*.jpg;*.png;*.bmp;*.tif'pecigy the extension of your image file
for i = 1:numel(files);
  filename = files(i).name;
  image = imread(filename);
  var=var+1;
  x=x+2;
  
  % apply processing to the loaded image
  % save the image to an array or back to the folder using 
end
if(var~=0)
fprintf('\nVar:%d\nX:%d\n',var,x); 
end
if(a)
    disp('True');
end
