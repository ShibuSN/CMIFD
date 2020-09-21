clc;
clear;
close all;
addpath('/home/shibuzzz/Downloads/vlfeat-0.9.21-bin/vlfeat-0.9.21/toolbox');
vl_setup; 
show=1;
color=[64,0,128];
TR_P=2;
TR_sim=15;
%Get File
[name, path]=uigetfile({'*.jpg;*.png;*.bmp;*.tif','Images (*.jpg,*.png,*.bmp,*.tif)'},'Select An Image');
if isequal(name,0)
    error('User selected Cancel');
end

filename=fullfile(path, name);
RGBimage=imread(filename);
tic;
[M,N,~]=size(RGBimage);
grayimage=rgb2gray(RGBimage);
Final_map1=false(M,N);
Final_map2=false(M,N);
Final_map3=false(M,N);

%% Segmentation
[CA1,CB1,CC1,CD1]=dwt2(RGBimage,'haar');
[CA2,CB2,CC2,CD2]=dwt2(CA1,'haar');
[CA3,CB3,CC3,CD3]=dwt2(CA2,'haar');
[CA4,CB4,CC4,CD4]=dwt2(CA3,'haar');
E_LF=sum(abs(CA4(:)));
E_HF=sum(abs(CB1(:)))+sum(abs(CB2(:)))+sum(abs(CB3(:)))+sum(abs(CB4(:)))+...
    sum(abs(CC1(:)))+sum(abs(CC2(:)))+sum(abs(CC3(:)))+sum(abs(CC4(:)))+...
    sum(abs(CD1(:)))+sum(abs(CD2(:)))+sum(abs(CD3(:)))+sum(abs(CD4(:)));
P_LF=E_LF/(E_LF+E_HF)*100;
S=(P_LF>50)*sqrt(0.02*M*N)+(P_LF<=50)*sqrt(0.01*M*N);
segments=vl_slic(im2single(RGBimage),S,1)+1;
num_segments=max(segments(:));
if show==1
    boundaries=imfilter(segments,fspecial('laplacian'))>0.1;
    boundaries=repmat(boundaries,[1,1,3]);
    boundaries=imdilate(boundaries,ones(5));
    imp = RGBimage;
    imp(boundaries) = 0 ;
    figure;
    imshow(imp);
end
%% SIFT
[Locations,Descriptors]= vl_sift(single(grayimage));
Locations=Locations([2,1,3,4],:)';
Descriptors=Descriptors';   %Tanspose of the matrix 
% [~,Descriptors,Locations]=sift(filename);
num_keypoint=size(Locations,1);
%Descriptor Normalization
Descriptors=double(Descriptors);
Descriptors=Descriptors./repmat(NormRow(Descriptors,2),1,128);
% Descriptors=Descriptors./vecnorm(Descriptors,2);
%% Matching
% Assign a segment number to each keypoint
keypoint_segment=zeros(num_keypoint,1);
for k=1:num_keypoint
    keypoint_segment(k)=segments(round(Locations(k,1)),round(Locations(k,2)));
end
%STEP1: Compare each 2 segments and create a matchList
CC=zeros(num_segments*(num_segments-1)/2,1);
MatchList=[];
num=0;
for i=1:num_segments
    idx1=find(keypoint_segment==i);
    if length(idx1)<2
        num=num+num_segments-i;
        continue;
    end
    D1=Descriptors(idx1,:);
    min1=min(pdist(D1));
    for j=i+1:num_segments
        num=num+1;
        idx2=find(keypoint_segment==j);
        if length(idx2)<2
            continue;
        end
        D2=Descriptors(idx2,:);
        min2=min(pdist(D2));
        min3=min(min1,min2);
        feat_dist=pdist2(D1,D2)/min3;%Distances between descriptors
        [k1,k2]=find(feat_dist<=1/TR_P);
        new_pairs=[idx1(k1),idx2(k2)]; %matched pairs
        MatchList=[MatchList;new_pairs];
        CC(num)=size(new_pairs,1);
    end
end
num_match=size(MatchList,1);
if num_match==0
    disp(' The Image Is Not Forged..!!');
    return;
end
%Show Match Keypoints
if (show==1)
    figure;
    imshow(RGBimage);
    hold on;
    title('Match Keypoints');
    loc1_temp=Locations(MatchList(:,1),1:2);
    loc2_temp=Locations(MatchList(:,2),1:2);
    locs_temp=[loc1_temp;loc2_temp];  
    plot(locs_temp(:,2),locs_temp(:,1),'mo');%Show points
    temp=[loc1_temp,loc2_temp,nan(size(loc2_temp,1),4)]';
    locs_temp=reshape(temp(:),4,[]);
    plot([locs_temp(2,:);locs_temp(4,:)],[locs_temp(1,:);locs_temp(3,:)],'b-');
    hold off;
    width =100;
    height =100;
    for i=1:num_match
    x1Center = loc1_temp(i,2);
    y1Center = loc1_temp(i,1);
    x1 = x1Center - width/2;
    y1 = y1Center - height/2;
    x2Center = loc2_temp(i,2);
    y2Center = loc2_temp(i,1);
    x2 = x2Center - width/2;
    y2 = y2Center - height/2;
    posi=rectangle('Position', [xLeft,yBottom,width,height],'EdgeColor','b','LineWidth',1);  
    RGBimage1=imcrop(RGBimage,[x1,y1,width,height]);
    RGBimage2=imcrop(RGBimage,[x2,y2,width,height]);
    % Get size of existing image A.
    [rowsA,colsA,numberOfColorChannelsA] = size(RGBimage1);
    % Get size of existing image B.
    [rowsB,colsB,numberOfColorChannelsB] = size(RGBimage2);
    if rowsB ~= rowsA || colsA ~= colsB
    % Size of B does not match A, so resize B to match A's size.
    RGBimage2 = imresize(RGBimage2, [rowsA colsA]);
    end
    RGBimageFinal = cat(2,RGBimage1,RGBimage2);
    figure;
    imshow(RGBimage1);
    imshow(RGBimage2);
    hold on;
    addpath('Zernike Moment');
    method='';%name of method for saving result
    % initilization
%     resultName=strcat('results\',method,filename);
    color=[0,0,255];
    tic;%record time for calculate elapsed time
    grayimage1=double(rgb2gray(RGBimageFinal)); % convert RGB to gray
%     grayimage2=double(rgb2gray(RGBimage2));
    
%     finalImg=cat(2,grayimage1,grayimage2); 
%     [M ,N]=size(grayimage);
    B=10;%Block Dimension =B x B
    search_th=50;%threshold for search length in matching
    distance_th=50;
%     num_blocks=(M-B+1)*(N-B+1);%number of blocks
    FVsize=12;%Feature Vector Length
    Similarity_threshold=0.1;%factor of similarity
    % get feature matrix (extract feature of all overlapped blocks)
    [FeatureMatrix,Locations]=getFeatureMatrix_parallel(grayimage1,B,FVsize);
%     [FeatureMatrix2,Locations2]=getFeatureMatrix_parallel(grayimage2,B,FVsize);
%     [FeatureMatrix,Locations]=getFeatureMatrix_parallel(grayimage2,B,FVsize);
    % Sorting
%     FeatureMatrix=cat(1,FeatureMatrix1,FeatureMatrix2);
    [FeatureMatrix,index]=sortrows(FeatureMatrix);
%     Locations=cat(1,Locations1,Locations2);
    Locations=Locations(index,:);
    % Matching (finding similar blocks)
    MatchList=getMatches_parallel(FeatureMatrix,Locations,Similarity_threshold,search_th,distance_th);
    num_matches=size(MatchList,1);
    disp(num_matches);
%     showMatches(RGBimage,MatchList,Locations,B,color,1,1,'FINAL IMG');
    end
    
    
%     roi = images.roi.Rectangle(gca,'Position',[26.331,77.041,80,100]);
%     h = images.roi.Circle(gca,'Center',[26.331 77.041],'Radius',50);
% width =100;
% height =100;
% xCenter = 253.321;
% yCenter = 135.485;
% xLeft = xCenter - width/2;
% yBottom = yCenter - height/2;
% posi=rectangle('Position', [xLeft,yBottom,width,height],'EdgeColor','b','LineWidth',1);  
% RGBimage1=imcrop(RGBimage,[xLeft,yBottom,width,height]);
% figure;
% imshow(RGBimage1);    
% hold on;
% grayimage=double(rgb2gray(RGBimage1)); % convert RGB to gray
end
%STEP2: Calculate the TR_B
CC_Sorted=sortrows(CC);
der_1st=imfilter(CC_Sorted,[-1;1],'replicate');
mean_der_1st=mean(der_1st);
der_2nd=imfilter(CC_Sorted,[1;-2;1],'replicate');
TR_B=min(der_2nd(der_2nd>mean_der_1st));
%STEP3: Removing blocks with low matched points
CC_Square=squareform(CC)>=TR_B;
%STEP4: removing matched points such that CC_Square==0
select=false(1,num_match);
for k=1:num_match
    s1=keypoint_segment(MatchList(k,1));
    s2=keypoint_segment(MatchList(k,2));
    select(k)=CC_Square(s1,s2);
end
new_MatchList=MatchList(select,:);
num_match=size(new_MatchList,1);
if num_match==0
    disp('low matched points....');
    return;
else
    disp('Image is forged...!!!');
end
%% Forgery Region Extraction
%STEP1: new segmentation and assign segments to each matched keypoint
superpixel_size=(max(M,N)>=1500)*10+10;
superpixel=vl_slic(im2single(RGBimage),superpixel_size,1)+1;
num_superpixel=max(superpixel(:));
% Assign a superpixel number to each matched keypoint
Match_superpixel=zeros(num_match,2);%SR
superpixel_select=false(num_superpixel,1);%Which superpixels have match and must processed
for k=1:num_match
    t1=new_MatchList(k,1);
    t2=new_MatchList(k,2);
    Match_superpixel(k,1)=superpixel(round(Locations(t1,1)),round(Locations(t1,2)));
    Match_superpixel(k,2)=superpixel(round(Locations(t2,1)),round(Locations(t2,2)));
    superpixel_select(Match_superpixel(k,1))=1;
    superpixel_select(Match_superpixel(k,2))=1;
end
superpixel_select_list=find(superpixel_select);%Convert logical vector to index
superpixel_neighbor=zeros(num_superpixel,8);
superpixel_center=zeros(num_superpixel,2);
superpixel_color=zeros(num_superpixel,1);
%STEP2: find neighbors of matched segments
%computing the center of each superpixel
for i=1:num_superpixel
    map=(superpixel==i);
    [y,x]=find(map);
    superpixel_center(i,1)=mean(y);
    superpixel_center(i,2)=mean(x);
end
%finding neighbors + create initial final_map
for i=1:length(superpixel_select_list)
    map=(superpixel==superpixel_select_list(i));
    Final_map1(map)=1;
    %Finding neighbor superpixels
    center=superpixel_center(superpixel_select_list(i),:);
    dil_map=imdilate(map,true(7));
    [y,x]=find(dil_map-map>0);
    neighbors=false(num_superpixel,1);
    for j=1:length(x)
        temp=superpixel(y(j),x(j));
        neighbors(temp)=1;
    end
    idx1=find(neighbors);
    neighbor_centers=superpixel_center(idx1,:);
    neighbor_angles=atan2(neighbor_centers(:,1)-center(1),neighbor_centers(:,2)-center(2))*180/pi;
    for j=1:8
        idx2=find(difference_angular(neighbor_angles,j*45-180)<=22);
        if isempty(idx2)
            continue;
        end
        superpixel_neighbor(superpixel_select_list(i),j)=idx1(idx2(1));
    end
end
%Compute mean color of found neighbor superpixels
temp=false(num_superpixel,1);
temp(superpixel_neighbor(superpixel_neighbor~=0))=1;
temp=find(temp);
for i=1:length(temp)
    map=repmat(superpixel==temp(i),[1,1,3]);
    superpixel_color(temp(i))=mean(RGBimage(map));
end
%Compare neighbors using its colors
Final_map2=Final_map1;
for i=1:num_match
    s1=Match_superpixel(i,1);
    s2=Match_superpixel(i,2);
    for j=1:8
        n1=superpixel_neighbor(s1,j);
        n2=superpixel_neighbor(s2,j);
        if n1==0||n2==0
            continue;
        end
        if abs(superpixel_color(n1)-superpixel_color(s1))<=TR_sim || abs(superpixel_color(n2)-superpixel_color(s2))<=TR_sim
            Final_map2(superpixel==n1|superpixel==n2)=1;
        end
    end
end
%STEP3: Morphological Operation
Final_map3=imclose(Final_map2,getCircleMask(2*superpixel_size));
toc
%% Create Result
% figure;imshow(map2RGB(Final_map1,RGBimage,color));title('Before Merging Neighbors');
% figure;imshow(map2RGB(Final_map2,RGBimage,color));title('After Merging Neighbors');
% figure;imshow(map2RGB(Final_map3,RGBimage,color));title('After Morphology');