classdef RewWarp
% Brief: 实现REW（Robust Elastic Warpin）图像变形算法，用于在存在视差
%   或视角差异时对移动图像进行弹性畸变校正以实现更自然的拼接结果。
% Details:
%   本类基于论文：
%     J. Li, Z. Wang, S. Lai, Y. Zhai and M. Zhang, "Parallax-Tolerant Image
%     Stitching Based on Robust Elastic Warping," IEEE Trans. on Multimedia,
%     vol.20, no.7, pp.1672-1687, Jul. 2018.
%   REW 在局部使用 Thin-Plate Spline (TPS) 拟合 residual (fixed - moving)
%   并结合到全局变换的平滑过渡，从而在保持局部几何一致性的同时
%   缓解视差带来的错位。
% 
%
% Syntax:
%   obj = RewWarp(movingPts, fixedPts)
%   obj = RewWarp(movingPts, fixedPts, lambda)
%
% Inputs:
%   movingPts  - Nx2 double, 源控制点坐标 (例如 img2 中的点)，每行 [x, y]
%   fixedPts   - Nx2 double, 目标控制点坐标 (仍然对应 img2 中的点)，每行 [x, y]
%   lambda     - (optional) double, TPS 正则化权重（默认 3e3），在构造函数中
%                内部乘以 8*pi 以匹配实现的数值尺度
%
% Outputs:
%   obj        - RewWarp 对象，包含 TPS 参数 (wx, wy, ax, bx, cx, ay, by, cy)
%                以及用于计算残差和对图像进行变形的方法
%
% Public methods (主要接口):
%   inverseMapping(obj, fixedPts)
%       对给定的目标坐标 fixedPts 计算 residual = fixed - moving 的位移，
%       返回 Nx2 的 [dx, dy]。
%
%   [warped,gx,hy,eta] = warpImage(obj, img2, meshImg2_U, meshImg2_V, meshPano_U, meshPano_V, offset_x, offset_y, overlapLimitsU, overlapLimitsV, K_smooth, interval_Mesh)
%       将移动图像 img2 根据 REW 变形并采样到 pano 网格上。返回变形后图像、
%       pano 网格上的残差 gx, hy 以及过渡权重 eta。方法内部在 overlap 区域
%       计算 TPS 残差并与全局变换平滑混合以获得最终采样坐标。
%
% Example:
%   见REW.mlx调用示例。
%
% Notes / 注意:
%   - 所有坐标均采用 MATLAB 像素坐标（1-based，格式 [x,y]）。
%   - 构造函数拟合的是 fixedPts - movingPts 的残差；之后通过权重分布进行
%     迭代剔除异常点（3-sigma 规则）。
%   - 参数 lambda 越大，拟合越平滑（更强的正则化）；根据点数与噪声可调节。
% 
% See also: warpImageTPS,warpImageTPS3,tpaps,fit

% Author:                          cuixingxing
% Email:                           cuixingxing150@gmail.com
% Created:                         20-Nov-2025 09:21:14
% Version history revision notes:
%                                  None
% Implementation In Matlab R2024a
% Copyright © 2025 xingxingcui.All Rights Reserved.
%
    
    properties(SetAccess=private)
        movingPts  (:,2) double % Nx2 矩阵 (源点, e.g., img2 的点)
        fixedPts   (:,2) double% Nx2 矩阵 (目标点, e.g., img2 对应点)
    end

    properties(SetAccess=private)
        wx           (:,1) double      % TPS 权重 x 方向 (n×1)
        wy           (:,1) double      % TPS 权重 y 方向 (n×1)
        ax           (1,1) double      % 仿射 a 系数 for x
        bx           (1,1) double      % 仿射 b 系数 for x
        cx           (1,1) double      % 常数项 for x
        ay           (1,1) double      % 仿射 a 系数 for y
        by           (1,1) double      % 仿射 b 系数 for y
        cy           (1,1) double      % 常数项 for y
        lambda       (1,1) double      % 正则化系数
         
        gxy          (:,2) double      % 对应控制点的残差,形如[x,y]两列，大小为nx2
        inlierIdx    (:,1) double      % 对应原始输入控制点（内点）的索引
    end
    
    methods(Access=public)
        function obj = RewWarp(movingPts, fixedPts, lambda)
            % Constructor：建立 TPS 模型
            arguments
                movingPts (:,2) double
                fixedPts (:,2) double
                lambda (1,1) double = 3e3; % weighting parameter to balance the fitting term and the smoothing term
            end
            obj.lambda = lambda* 8*pi;
            
            % 构建 TPS 系数（w, a, b, c）
            % 这里我们是对fixedPts-movingPts 的残差拟合 
            n = size(fixedPts,1);
            X = fixedPts(:,1);
            Y = fixedPts(:,2);

            % 计算 K 矩阵
            dist2 = (X - X').^2 + (Y - Y').^2;
            % 避免 log(0)
            mask = dist2 == 0;
            dist2(mask) = 1;
            K = 0.5 * dist2 .* log(dist2);
            K(mask) = 0;

            % 正则化
            K = K + obj.lambda * eye(n);
            
            % P 矩阵 (affine 部分)
            P = [X, Y,ones(n,1)];
            
            % 目标 residual
            G =  fixedPts - movingPts; % Nx2
            
            % 构建线性系统
            A = [K, P; P', zeros(3,3)];
            Gxy = [G; zeros(3,2)];
            w = A\Gxy;
            
            wxa = w(:,1);
            mwa = w(:,2);
            
            obj.wx = wxa(1:n);
            obj.ax = wxa(n+1);
            obj.bx = wxa(n+2);
            obj.cx = wxa(n+3);
            
            obj.wy = mwa(1:n);
            obj.ay = mwa(n+1);
            obj.by = mwa(n+2);
            obj.cy = mwa(n+3);

            % remove outliers based on the distribution of weights，即3sigma准则排除
            outlier = abs(obj.wx)>3*std(obj.wx) | abs(obj.wy)>3*std(obj.wy);

            inlier_idx = 1:n;
            for kiter = 1:10
                if sum(outlier) < 0.0027*n
                    break;
                end
                ok = ~outlier;
                inlier_idx = inlier_idx(ok);
                A = A([ok;true(3,1)],[ok;true(3,1)]);
                Gxy = Gxy([ok;true(3,1)],:);
                W = A\Gxy;
                n = length(inlier_idx);

                obj.wx = W(1:n,1);
                obj.wy = W(1:n,2);

                obj.ax = W(n+1,1);
                obj.bx = W(n+2,1);
                obj.cx = W(n+3,1);

                obj.ay = W(n+1,2);
                obj.by = W(n+2,2);
                obj.cy = W(n+3,2);

                outlier = abs(obj.wx)>3*std(obj.wx) | abs(obj.wy)>3*std(obj.wy);  
            end

            obj.movingPts = movingPts(inlier_idx,:);
            obj.fixedPts = fixedPts(inlier_idx,:);
            obj.gxy = G(inlier_idx,:);
            obj.inlierIdx = inlier_idx;
        end
        

        
        function [warped,gx,hy,eta] = warpImage(obj, img2, meshImg2_U, meshImg2_V, meshPano_U,meshPano_V,offset_x,offset_y,overlapLimitsU,overlapLimitsV,K_smooth,interval_Mesh)
            % 对图像 img (moving 图像) 使用 REW warp (TPS + 混合)
            arguments
                obj
                img2 {mustBeNumeric}
                meshImg2_U (:,:) double% 经global homography/similarity转换后对应的图像2网格横坐标
                meshImg2_V (:,:) double% 经global homography/similarity转换后对应的图像2网格纵坐标
                meshPano_U (:,:) double% 经global homography/similarity转换后对应的全景图像网格横坐标
                meshPano_V (:,:) double% 经global homography/similarity转换后对应的全景图像网格纵坐标
                offset_x (1,1) double % meshImg2_U(1)对应meshPano_U中的横向偏置，即x横坐标
                offset_y (1,1) double % meshImg2_V(1)对应meshPano_V中的纵向偏置，即y纵坐标
                overlapLimitsU (1,2) double% 经global homography/similarity转换后对应的overlap图像网格横坐标范围,形如[umin,umax]
                overlapLimitsV (1,2) double% 经global homography/similarity转换后对应的overlap图像网格纵坐标范围,形如[vmin,vmax]
                K_smooth (1,1) double = 5 % 光滑程度阈值，数值越大，过渡越光滑
                interval_Mesh (1,1) double = 10 % 插值网格采样点间隔，数值越大，采样点越稀疏，计算效率越高
            end

            %% 全局所有网格点残差TPS计算
            [meshH,meshW] = size(meshImg2_U);
            meshU = meshImg2_U(1:interval_Mesh:end,1:interval_Mesh:end);
            meshV = meshImg2_V(1:interval_Mesh:end,1:interval_Mesh:end);

            meshFixedPts = [meshU(:),meshV(:)];
            residualPts = inverseMapping(obj,meshFixedPts);

            % 恢复采样网格点形状
            u_residual = reshape(residualPts(:,1),size(meshU));
            v_residual = reshape(residualPts(:,2),size(meshV));
            gx_sub = imresize(u_residual, [meshH,meshW]);
            hy_sub = imresize(v_residual, [meshH,meshW]);

            % 全景图残差gx,gy分配
            [panoHeight,panoWidth] = size(meshPano_U);
            gx = zeros(panoHeight,panoWidth);
            hy = zeros(panoHeight,panoWidth);

            x_span =  offset_x:offset_x+meshW-1; 
            y_span = offset_y:offset_y+meshH-1;

            % 避免越界
            x_span(x_span<1) = 1;
            x_span(x_span>panoWidth) = panoWidth;
            y_span(y_span<1) = 1;
            y_span(y_span>panoHeight) = panoHeight;

            gx(y_span,x_span) = gx_sub;
            hy(y_span,x_span) = hy_sub;

            %% smooth tansition to global transform
            eta_d0 = 0; % lower boundary for smooth transition area
            eta_d1 = K_smooth * max(abs(obj.gxy),[],"all"); % higher boundary for smooth transition area

            % [xLimitsOutInImg2,yLimitsOutInImg2] = outputLimits(tform,[1,size(img1,2)],[1,size(img1,1)]);
            sub_u0_ = max(1,overlapLimitsU(1));
            sub_u1_ = min(size(img2,2),overlapLimitsU(2));
            sub_v0_ = max(1,overlapLimitsV(1));
            sub_v1_ = min(size(img2,1),overlapLimitsV(2));

            sub_u0_ = sub_u0_ + min(obj.gxy(:,1)); % 这里作用不大，也可以去掉
            sub_u1_ = sub_u1_ + max(obj.gxy(:,1));% 这里作用不大，也可以去掉
            sub_v0_ = sub_v0_ + min(obj.gxy(:,2));% 这里作用不大，也可以去掉
            sub_v1_ = sub_v1_ + max(obj.gxy(:,2));% 这里作用不大，也可以去掉

            dist_horizontal = max(sub_u0_-meshPano_U, meshPano_U-sub_u1_);
            dist_vertical = max(sub_v0_-meshPano_V, meshPano_V-sub_v1_);
            dist_sub = max(dist_horizontal, dist_vertical);
            dist_sub = max(0, dist_sub);
            eta = (eta_d1 - dist_sub) ./ (eta_d1 - eta_d0); % gxy的加权系数，非overlapp区域逐渐衰减为零，大小同全景图
            eta(dist_sub < eta_d0) = 1;
            eta(dist_sub > eta_d1) = 0;
            
            gx = gx .* eta;
            hy = hy .* eta;

            % 计算最终采样坐标
            meshPano_U = meshPano_U - gx;
            meshPano_V = meshPano_V - hy;
            
            % 插值
            warped = images.internal.interp2d(img2,meshPano_U,meshPano_V,"linear",0, false);
        end
    end

    methods(Access=private)
        function residualPts = inverseMapping(obj,fixedPts)
            % fixedPts为矫正后的图像坐标，返回fixedPts-movingPts的残差坐标
            arguments
                obj 
                fixedPts (:,2) double 
            end
            u_im_ = fixedPts(:,1);
            v_im_ = fixedPts(:,2);

            % 构造目标网格到控制点的径向基
            npix = numel(u_im_);
            n = size(obj.fixedPts,1);
            rbf = zeros(npix,n);
            for kf = 1:n
                dist2 = (u_im_ - obj.fixedPts(kf,1)).^2 + (v_im_ - obj.fixedPts(kf,2)).^2;
                rbf(:,kf) = 0.5 * dist2 .* log(dist2);
            end

            % 目标网格点反向映射到源图像残差浮点坐标
            gx_sub = rbf * obj.wx +  obj.ax.*u_im_+obj.bx.*v_im_+obj.cx;
            hy_sub = rbf * obj.wy +  obj.ay.*u_im_+obj.by.*v_im_+obj.cy;

            residualPts = [gx_sub,hy_sub]; % 返回[npix,2]大小的残差坐标
        end
    end
end

