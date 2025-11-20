function [outImg, outR] = warpImageTPS2(img, movingPoints, fixedPoints,options)
% Brief: 基于图像上匹配的控制点使用 Thin Plate Spline 对图像做非刚性形变
% Details:
%   根据位于`img`中的二维点集从`movingPoints`移动到匹配对应的二维点集`fixedPoints`
% 关系，对图像`img`进行TPS非刚性变形，返回输出图像`outImg`（附带空间参考对象outR）
%
% Syntax:  
%   outImg = warpImageTPS2(img, pts1, pts2)
%   [outImg, outR] = warpImageTPS2(img, pts1, pts2, options)
% 
% Inputs:
%   img             - 输入图像，A
%   movingPoints    - 源控制点 Nx2,即输入图像img上需要移动的像素点坐标
%   fixedPoints     - 目标控制点 Nx2,即输入图像img上对应移动到的像素点坐标
%   options         - 可选参数结构体，字段如下
%       .refImg     - 图像空间参考对象，imref2d
%       .outputView - 输出图像空间参考对象，imref2d
%       .lambda     - 正则化参数（默认0）
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
%    outImg2 = warpImageTPS3(distortImg, distortPts, correctPts,outputView=imref2d(size(img)));
%    figure;
%    imshowpair(distortImg,outImg2,"montage")
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
    movingPoints (:,2) double
    fixedPoints (:,2) double
    options.refImg (1,1) imref2d = imref2d(size(img))
    options.outputView (1,1) imref2d = imref2d(size(img))
    options.lambda double = 0
end

refImg = options.refImg;
outputView = options.outputView;
lambda = options.lambda;

% 像素点转为世界坐标点
[fixedWorldPtsX, fixedWorldPtsY] = intrinsicToWorld(refImg, fixedPoints(:,1), fixedPoints(:,2));
fixedWorldPts = [fixedWorldPtsX,fixedWorldPtsY];
[movingPtsX, movingPtsY] = intrinsicToWorld(refImg, movingPoints(:,1), movingPoints(:,2));
movingWorldPts = [movingPtsX,movingPtsY];

% 1. 构建TPS线性系统: 输入为目标点(pts1)，输出为源点(pts2)
n = size(fixedWorldPts,1);
X = fixedWorldPts(:,1);
Y = fixedWorldPts(:,2);

% 构造核矩阵
K = zeros(n, n);
for i = 1:n
    for j = 1:n
        r2 = (X(i)-X(j))^2 + (Y(i)-Y(j))^2;
        if r2==0
            K(i,j) = 0;
        else
            K(i,j) = r2 * log(r2) * 0.5;
        end
    end
end
% 正则化
K = K + lambda*eye(n);

P = [ones(n,1), X, Y];

% 右端项(pts2，目标点→源点插值)
v1 = movingWorldPts(:,1); % x
v2 = movingWorldPts(:,2); % y

% 拼合线性系统
A = [K, P; 
    P', zeros(3,3)];
vx = [v1; zeros(3,1)];
vy = [v2; zeros(3,1)];

% 求解参数
params_x = A \ vx;
params_y = A \ vy;

w_x = params_x(1:n);
a_x = params_x(n+1:end);
w_y = params_y(1:n);
a_y = params_y(n+1:end);

% 2. 对目标输出网格做反向映射
[xq, yq] = meshgrid(1:outputView.ImageSize(2), 1:outputView.ImageSize(1));
[xWorld, yWorld] = intrinsicToWorld(outputView, xq, yq);
npix = numel(xWorld);

Xq = xWorld(:);
Yq = yWorld(:);

% 构造目标网格到控制点的径向基
U = zeros(npix, n);
for i = 1:n
    r2 = (Xq - X(i)).^2 + (Yq - Y(i)).^2;
    U(:,i) = 0.5 * r2 .* log(max(r2,eps));
end

% 目标网格点反向映射到源图像浮点坐标
uSource = U * w_x + a_x(1) + a_x(2)*Xq + a_x(3)*Yq;
vSource = U * w_y + a_y(1) + a_y(2)*Xq + a_y(3)*Yq;

% 转为输入图像内部坐标
[uSourceI, vSourceI] = worldToIntrinsic(refImg, uSource, vSource);
uSourceI = reshape(uSourceI, size(xq));
vSourceI = reshape(vSourceI, size(yq));

% 3. 双线性插值采样
outImg = imageInterp(img,uSourceI, vSourceI);
outR = outputView;
end