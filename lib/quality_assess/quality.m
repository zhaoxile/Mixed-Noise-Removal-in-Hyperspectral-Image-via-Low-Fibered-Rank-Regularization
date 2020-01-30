function [psnr] = quality(imagery1, imagery2)
%==========================================================================
% Evaluates the quality assessment indices for two tensors.
%
% Syntax:
%   [psnr, ssim, fsim] = quality(imagery1, imagery2)
%
% Input:
%   imagery1 - the reference tensor
%   imagery2 - the target tensor

% NOTE: the tensor is a M*N*K array and DYNAMIC RANGE [0, 255]. 

% Output:
%   psnr - Peak Signal-to-Noise Ratio
%   ssim - Structure SIMilarity
%   fsim - Feature SIMilarity

% See also StructureSIM, FeatureSIM
%
% by Yi Peng
% update by Yu-Bang Zheng 11/19/2018
%==========================================================================
Nway = size(imagery1);
psnr = zeros(Nway(3),1);
for i = 1:Nway(3)
    psnr(i) = psnr_index(imagery1(:, :, i), imagery2(:, :, i));
end


