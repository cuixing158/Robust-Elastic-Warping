function outImage = imageInterp(inImage,mapX,mapY,options)
% Brief: Interpolate image based on mapping coordinates
%
% Syntax:
%     outImage = imageInterp(inImage,mapX,mapY)
%     outImage = imageInterp(inImage,mapX,mapY,name=value)
%
% Inputs:
%   inImage   - 输入图像，大小为 [h, w, c]，支持灰度或彩色
%   mapX      - 目标像素的 X 坐标映射，大小为 [oH, oW]，类型为 double/single
%   mapY      - 目标像素的 Y 坐标映射，大小为 [oH, oW]，类型为 double/single
%   options   - 结构体或参数集，包含以下字段:
%       FillValues   - 边界外填充值（double，默认0），用于 BORDER_CONSTANT
%       BorderMode   - 边界模式（string,默认BORDER_CONSTANT），可选:
%                         "BORDER_CONSTANT"   : 超出边界用 FillValues 填充
%                         "BORDER_REPLICATE"  : 超出边界用最近边界像素值填充
%                         "BORDER_REFLECT"    : 超出边界做镜像反射
%       SmoothEdges  - 是否对边缘做平滑/反锯齿处理（logical，默认 false）
%       InterpMethod - 像素插值方法（string,默认"linear"），可选:
%                         "linear"  :线性插值
%                         "nearest" :最近邻插值
%                         "cubic"   :cubic插值
%
% Outputs:
%   outImage  - 变形后的输出图像，大小与mapX/mapY一致，类型与输入inImage一致
%
% See also: interp2, imwarp

% Author:                          cuixingxing
% Email:                           cuixingxing150@gmail.com
% Created:                         03-Jul-2025 20:32:13
% Version history revision notes:
%                                  17-Jul-2025 add interpolation methods 
% Implementation In Matlab R2025a
% Copyright © 2025 TheMatrix.All Rights Reserved.
%
arguments
    inImage {mustBeNumeric,mustBeFinite} % 输入图像
    mapX (:,:) double {mustBeNumeric,mustBeFinite} % 输入图像的 X 坐标映射，大小为 [oH, oW]，类型为 double/single
    mapY (:,:) double {mustBeNumeric,mustBeFinite} % 输入图像的 Y 坐标映射，大小为 [oH, oW]，类型为 double/single
    options.FillValues (1,1) double {mustBeInRange(options.FillValues,0,255)}=0
    options.BorderMode (1,1) string {mustBeMember(options.BorderMode,["BORDER_CONSTANT","BORDER_REPLICATE","BORDER_REFLECT"])} = "BORDER_CONSTANT"
    options.SmoothEdges logical = false; % 是否对边缘反锯齿/光滑处理
    options.InterpMethod (1,1) string {mustBeMember(options.InterpMethod,["linear","nearest","cubic"])} = "linear"
end
[h,w,~] = size(inImage);

switch options.BorderMode
    case "BORDER_CONSTANT"
        % 使用固定值填充边界
        outImage = images.internal.interp2d(inImage,mapX,mapY,options.InterpMethod,options.FillValues,options.SmoothEdges);

    case "BORDER_REPLICATE"
        % 限制坐标到图像范围内，实现边界复制
        mapX_clamped = max(1, min(w, mapX));
        mapY_clamped = max(1, min(h, mapY));
        % out = interp2(1:w, 1:h, inImage, mapX_clamped, mapY_clamped, 'linear');
        outImage = images.internal.interp2d(inImage,mapX_clamped,mapY_clamped,options.InterpMethod,options.FillValues,options.SmoothEdges);

    case "BORDER_REFLECT"
        % 反射坐标处理
        mapX_reflect = reflectCoordinates(mapX, w);
        mapY_reflect = reflectCoordinates(mapY, h);
        % out = interp2(1:w, 1:h, inImage, mapX_reflect, mapY_reflect, 'linear');
        outImage = images.internal.interp2d(inImage,mapX_reflect,mapY_reflect,options.InterpMethod,options.FillValues,options.SmoothEdges);
    otherwise
        error('Invalid border mode: %s', borderMode);
end
end

% 反射坐标处理
function coord = reflectCoordinates(coord, maxVal)
% 将坐标反射到1到maxVal的范围内
coord = 2*maxVal - coord;
coord = mod(coord - 1, 2*(maxVal - 1)) + 1;
coord = maxVal - abs(maxVal - coord);
coord = max(1, min(maxVal, coord));
end