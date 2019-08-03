function [nfm, phs, mask_sharp] = m1_bb_vSHARP(nfm_mac,msk,Kernel_Sizes)
    % Fast implementation of SHARP filtering with variable kernel size
    % Script originally written by B. Bilgic
    % For academic use only, please refer to the appropriate publications: 
    % https://www.martinos.org/~berkin/index.html
    %
    % Modified and optimized by: 
    % Job G. Bouwman
    % 25-11-2014
    % jgbouwman@hotmail.com
    % 
    % Copyright (c) 2016, Job
    % All rights reserved.
    % 
    % Redistribution and use in source and binary forms, with or without
    % modification, are permitted provided that the following conditions are
    % met:
    % 
    %     * Redistributions of source code must retain the above copyright
    %       notice, this list of conditions and the following disclaimer.
    %     * Redistributions in binary form must reproduce the above copyright
    %       notice, this list of conditions and the following disclaimer in
    %       the documentation and/or other materials provided with the distribution
    % 
    % THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
    % AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
    % IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
    % ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT OWNER OR CONTRIBUTORS BE
    % LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR
    % CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF
    % SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS
    % INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN
    % CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE)
    % ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE
    % POSSIBILITY OF SUCH DAMAGE.

    %% Embedding the matrices in an efficient Field of View:
        minimalBorder = 2+max(Kernel_Sizes);
        Nold = size(nfm_mac);
    
    % The original indices in which the ROI is located: 
        x1 = find(squeeze(sum(sum(msk, 1),3)));
        y1 = find(squeeze(sum(sum(msk, 2),3)));
        z1 = find(squeeze(sum(sum(msk, 1),2)));

    % The size of the new matrix:
        NyxzNew = 4*ceil(([length(y1),length(x1),length(z1)] + minimalBorder)/4);

    % The new indices in which the ROI will be located: 
        x2 = x1 + (-x1(1) + round((NyxzNew(2) - length(x1))/2));
        y2 = y1 + (-y1(1) + round((NyxzNew(1) - length(y1))/2));
        z2 = z1 + (-z1(1) + round((NyxzNew(3) - length(z1))/2));
    % Embedding the new matrices:
        nfm_macNew = zeros(NyxzNew); 
        nfm_macNew(y2,x2,z2) = nfm_mac(y1, x1, z1); 
        nfm_mac = nfm_macNew; clear nfm_macNew;
        
        maskNew = zeros(NyxzNew); 
        maskNew(y2,x2,z2) = msk(y1, x1, z1); 
        msk = maskNew; clear maskNew;
        N = NyxzNew;
        
    %% Hacking the imaginary channel, the length of kern
        if mod(length(Kernel_Sizes), 2)
            Kernel_Sizes = [(2*Kernel_Sizes(1)-Kernel_Sizes(2)) Kernel_Sizes];
        end

    %% initialize
        threshold = .05;
        FT_nfm_mac = fftn(nfm_mac); clear nfm_mac;
        phs = 0; 
        mask_sharp_0 = 0;
        
    %% Forward SHARP filtering with variable kernel size:
    for m = 1:2:length(Kernel_Sizes) % (performed per pair)
        for n = 0:1
            Kernel_Size = Kernel_Sizes(m+n);
            ksize = [Kernel_Size, Kernel_Size, Kernel_Size];                % Sharp kernel size
%             khsize = (ksize-1)/2;
            khsize = round((ksize)/2);
            [a,b,c] = meshgrid(-khsize(2):khsize(2), -khsize(1):khsize(1), -khsize(3):khsize(3));
            kernel = (a.^2 / khsize(1)^2 + b.^2 / khsize(2)^2 + c.^2 / khsize(3)^2 ) <= 1;
            kernel = -kernel / sum(kernel(:));
            kernel(khsize(1)+1,khsize(2)+1,khsize(3)+1) = 1 + kernel(khsize(1)+1,khsize(2)+1,khsize(3)+1);

            if n == 0
                Kernel = zeros(N);
                Kernel( 1+N(1)/2 - khsize(1) : 1+N(1)/2 + khsize(1), 1+N(2)/2 - khsize(2) : 1+N(2)/2 + khsize(2), 1+N(3)/2 - khsize(3) : 1+N(3)/2 + khsize(3) ) = kernel;
                ksize_1 = ksize(1);
            else
                Kernel( 1+N(1)/2 - khsize(1) : 1+N(1)/2 + khsize(1), 1+N(2)/2 - khsize(2) : 1+N(2)/2 + khsize(2), 1+N(3)/2 - khsize(3) : 1+N(3)/2 + khsize(3) ) = ...
                Kernel( 1+N(1)/2 - khsize(1) : 1+N(1)/2 + khsize(1), 1+N(2)/2 - khsize(2) : 1+N(2)/2 + khsize(2), 1+N(3)/2 - khsize(3) : 1+N(3)/2 + khsize(3) ) + ...
                    + (1i)*kernel;
                ksize_2 = ksize(1);
            end
        end
        
        % (Here we 'hack' Matlab's FFT, using the imaginary channel as well)
        del_sharp_DUAL = fftn(ifftshift(Kernel)); 
            
        % creating the erosion masks for the two kernels: 
        erode_size = ksize_1 + 1;
        mask_sharp = imerode(msk       , strel('line', erode_size, 0));
        mask_sharp = imerode(mask_sharp, strel('line', erode_size, 90));
        mask_sharp = permute(mask_sharp, [1,3,2]);
        mask_sharp = imerode(mask_sharp, strel('line', erode_size, 0));
        mask_sharp_1 = permute(mask_sharp, [1,3,2]);

        erode_size = ksize_2 + 1;
        mask_sharp = imerode(msk       , strel('line', erode_size, 0));
        mask_sharp = imerode(mask_sharp, strel('line', erode_size, 90));
        mask_sharp = permute(mask_sharp, [1,3,2]);
        mask_sharp = imerode(mask_sharp, strel('line', erode_size, 0));
        mask_sharp_2 = permute(mask_sharp, [1,3,2]);
       
        % apply Sharp to normalized field shift:
        phs_dual    = ifftn(FT_nfm_mac.*del_sharp_DUAL);

        % cumulative sharp_filtered field
        phs = phs + (mask_sharp_1 - mask_sharp_0).*real(phs_dual) + ...
                    (mask_sharp_2 - mask_sharp_1).*imag(phs_dual);
                
        % for the next pair of masks:
        mask_sharp_0     = mask_sharp_2;
        
        if m == 1
            delsharp_inv = 1./real(del_sharp_DUAL);
            delsharp_inv(real(del_sharp_DUAL)<threshold) = 0;
        end
    end
    %% Inverse SHARP filtering:
    nfm = real( ifftn(fftn(phs).*delsharp_inv) .* mask_sharp_0); 
    
    %% reembedding into original matrix size:
    nfm_2 = zeros(Nold);
    nfm_2(y1, x1, z1) = nfm(y2, x2, z2);
    nfm = nfm_2;
    
    phs_2 = zeros(Nold);
    phs_2(y1, x1, z1) = phs(y2, x2, z2);
    phs = phs_2;
  
    mask_sharp = zeros(Nold);
    mask_sharp(y1, x1, z1) = mask_sharp_0(y2, x2, z2);
end
