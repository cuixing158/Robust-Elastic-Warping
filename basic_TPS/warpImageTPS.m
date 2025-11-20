function [outImg, outR] = warpImageTPS(img, movingPoints, fixedPoints,options)
% Brief: 基于图像上匹配的控制点使用 Thin Plate Spline 对图像做非刚性形变
% Details:
%   根据位于`img`中的二维点集从`movingPoints`移动到匹配对应的二维点集`fixedPoints`
% 关系，对图像`img`进行TPS非刚性变形，返回输出图像`outImg`（附带空间参考对象outR）
%
% Syntax:  
%   outImg = warpImageTPS(img, movingPoints, fixedPoints)
%   [outImg, outR] = warpImageTPS(img, movingPoints, fixedPoints, options)
% 
% Inputs:
%   img             - 输入图像，A
%   movingPoints    - 源控制点 Nx2,即输入图像img上需要移动的像素点坐标
%   fixedPoints     - 目标控制点 Nx2,即输入图像img上对应移动到的像素点坐标
%   options         - 可选参数结构体，字段如下
%       .refImg     - 图像空间参考对象，imref2d
%       .outputView - 输出图像空间参考对象，imref2d
%       .lambda     - 正则化参数（默认0）
%       .K_smooth   - 平滑边界系数（默认1）,[0,1]之间，0代表光滑最小二乘拟合，1代表插值
% 
% Outputs:
%   outImg     - 形变后图像
%   outR       - 输出图像空间参考对象（imref2d）
% 
% Example: 
%    %% step1,制作一副棋盘网格畸变图像
%    squareWidth = 50;% in pixels
%    rows = 3;
%    cols = 4;
%    img = checkerboard(squareWidth,rows,cols);
%    figure(1);imshow(img)
%    
%    % make a distortion image
%    height = squareWidth*2*rows;
%    width = squareWidth*2*cols;
%    fx = (width+height)/2;
%    fy = fx;
%    cx = width/2;
%    cy = height/2;
%    k = [fx 0 cx;
%        0 fy cy;
%        0 0 1];
%    radialDistortion = [-0.3361,15.8921,28.22]; 
%    cameraParams = cameraParameters("K",k,"RadialDistortion",radialDistortion);
%    distortImg = undistortImage(img,cameraParams,"cubic",OutputView="full",FillValues=0.5);
%    figure(2);imshow(distortImg)
%    
%    %% step2,畸变点和对应的非畸变点
%    % [srcPtX,srcPtY] = ginput(20);
%    distortPts = [18.173       31.393
%           64.771       12.408
%           138.98       3.0888
%           213.88       12.408
%           260.14       30.012
%           3.3303        85.93
%           47.857         77.3
%           139.67       72.123
%           231.14       76.265
%           275.32        85.93
%           3.6754       155.65
%           47.857       167.04
%           138.98       171.53
%           229.76       168.08
%           273.94       156.34
%           18.863       212.61
%           65.461       232.28
%           139.67       241.95
%           213.19        230.9
%           259.79       211.92];
%    
%    [x,y] = meshgrid(1:2*squareWidth:width+1,1:2*squareWidth:height+1);
%    x = x'; % 对应上门行优先排列
%    y = y';
%    correctPts = [x(:),y(:)];
%   
%    %% step3,TPS矫正畸变图像
%    correctedImg = warpImageTPS2(distortImg,distortPts, correctPts, outputView=imref2d(size(img)));
%    figure;
%    imshowpair(img,correctedImg)
% 
% See also: warpImageTPS,tpaps,fit

% Author:                          cuixingxing
% Email:                           cuixingxing150@gmail.com
% Created:                         29-Jun-2025 09:48:54
% Version history revision notes:
%                                  None
% Implementation In Matlab R2025a
% Copyright © 2025 TheMatrix.All Rights Reserved.
%
arguments
    img {mustBeNumeric}
    movingPoints (:,2) double % 位于输出outImg上的匹配二维像素点坐标
    fixedPoints (:,2) double % 位于输入图像img上匹配的二维像素点坐标
    options.refImg (1,1) imref2d = imref2d(size(img))          % 默认输入img对应的refImg
    options.outputView (1,1) imref2d = imref2d(size(img))
    options.K_smooth {mustBeInRange(options.K_smooth,0,1)}= 1 % 默认插值，严格经过插值点
end

refImg = options.refImg;
outputView = options.outputView;
K_smooth= options.K_smooth;

% 像素点转为世界坐标点
[fixedWorldPtsX, fixedWorldPtsY] = intrinsicToWorld(refImg, fixedPoints(:,1), fixedPoints(:,2));
fixedWorldPts = [fixedWorldPtsX,fixedWorldPtsY];
[movingPtsX, movingPtsY] = intrinsicToWorld(refImg, movingPoints(:,1), movingPoints(:,2));
movingWorldPts = [movingPtsX,movingPtsY];

% 构建 TPS 插值函数（从 pts1 → pts2）,inverse mapping 
fxy = tpaps(fixedWorldPts',movingWorldPts',K_smooth);% K_smooth=1,代表只最小化误差 → 精确插值;0代表只最小化弯曲能量 → 纯光滑函数，不考虑数据点位置

% 或者使用下述语句做TPS
% fx = fit(fixedWorldPts, movingWorldPts(:,1), 'thinplateinterp');
% fy = fit(fixedWorldPts, movingWorldPts(:,2), 'thinplateinterp');

% 构建输出网格（目标空间）
[xq, yq] = meshgrid(1:outputView.ImageSize(2), 1:outputView.ImageSize(1));
[xWorld, yWorld] = intrinsicToWorld(outputView, xq, yq);

% 反向映射目标空间点 → 源图像点
xySource = fnval(fxy,[xWorld(:)'; yWorld(:)']);
xSource = xySource(1,:);
ySource = xySource(2,:);

% 或者使用下述语句，对应上述fit函数对象
% xSource = fx([xWorld(:), yWorld(:)]);
% ySource = fy([xWorld(:), yWorld(:)]);

% 转换为源图像内部坐标（intrinsic）
[uSource, vSource] = worldToIntrinsic(refImg, xSource, ySource);
uSource = reshape(uSource, size(xq));
vSource = reshape(vSource, size(yq));

% 插值采样
outImg = imageInterp(img,uSource,vSource);
% outImg = zeros([outputView.ImageSize size(img,3)], class(img));
% for c = 1:size(img,3)
%     outImg(:,:,c) = interp2(double(img(:,:,c)), uSource, vSource, 'linear', 0);
% end

outR = outputView;
end
