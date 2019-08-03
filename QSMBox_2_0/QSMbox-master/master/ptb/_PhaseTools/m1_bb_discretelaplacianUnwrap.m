function [phaseUnwrapped phaseLaplacianFiltered] = m1_bb_discretelaplacianUnwrap(phaseWrapped, msk)
    % Laplacian unwrapping, with the Discrete Laplacian (DL) 
    %
    % Modification of equation 5 in Schofield et al, 
    % "Fast phase unwrapping algorithms for interferometric 
    % applications", OPTICS LETTERS, VOL. 28, No 14, 2003
    %
    % Written by: 
    % Job G. Bouwman
    % 25-11-2014
    % jgbouwman@hotmail.com
    %
    % First change: the forward discrete Laplacian is implemented instead
    % of the continuous Laplacian. (see: Bakker CJG, “Selective depiction 
    % of susceptibility transitions using Laplace-filtered phase maps.” 
    % Volume 30, Issue 5, June 2012, Pages 601–609
    %
    % Second change: correspondingly, the inverse Laplacian (deconvolution)
    % is also implemented with the discrete Laplacian
    %
    % Third change: the deconvolution part (fft's on real data, 
    % is performed by ffts on complex data.
    %
    % Usage: 
    % 'phaseWrapped' should be in interval (-pi, pi]   
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
    
    DISCRETE_LAPLACIAN = true;
    
    %% Embedding the matrix in a smaller FOV: 
        minimalBorder = 4;
        Nold = size(phaseWrapped);

        % The original indices in which the ROI is located: 
            x1 = find(squeeze(sum(sum(msk, 1),3)));
            y1 = find(squeeze(sum(sum(msk, 2),3)));
            z1 = find(squeeze(sum(sum(msk, 1),2)));

        % The size of the new matrix:
            NyxzNew = 4*ceil([length(y1),length(x1),length(z1)]/4) + 2*minimalBorder;

        % The new indices in which the ROI will be located: 
            x2 = x1 + (-x1(1) + round((NyxzNew(2) - length(x1))/2));
            y2 = y1 + (-y1(1) + round((NyxzNew(1) - length(y1))/2));
            z2 = z1 + (-z1(1) + round((NyxzNew(3) - length(z1))/2));
        % Embedding the new matrix:
            phaseWrappedNew = zeros(NyxzNew); 
            phaseWrappedNew(y2,x2,z2) = phaseWrapped(y1, x1, z1); 
            phaseWrapped = phaseWrappedNew; clear phaseWrappedNew;
    
    
    %% Preparation:
        % making the phaseWrapped data complex:
        P = exp(1i*phaseWrapped); clear phaseWrapped;
        N = size(P); 
        dkx = 1/N(2); 
        dky = 1/N(1); 
        dkz = 1/N(3); 
        
    %% Prepare the inverse (and forward) Laplacian filter:
    if DISCRETE_LAPLACIAN
        lap_x = 2 - 2*cos((0:(N(2)-1))*dkx*2*pi);
        lap_y = 2 - 2*cos((0:(N(1)-1))*dky*2*pi);
        lap_z = 2 - 2*cos((0:(N(3)-1))*dkz*2*pi);       
    else
        lap_x = (3/2)*(16)*ifftshift(((-N(2)/2:N(2)/2-1)*dkx).^2); 
        lap_y = (3/2)*(16)*ifftshift(((-N(1)/2:N(1)/2-1)*dky).^2); 
        lap_z = (3/2)*(16)*ifftshift(((-N(3)/2:N(3)/2-1)*dkz).^2);
    end
    [lap_x,lap_y,lap_z] = meshgrid(lap_x,lap_y,lap_z);
    del_op = (lap_x + lap_y + lap_z)/7; clear lap_x lap_y lap_z;
    del_inv = 1./del_op;
    del_inv(1,1,1) = 0;    
        
    
    %% Forward Laplacian phase filtering:
    if DISCRETE_LAPLACIAN
        % Start with first order derivative (complex division)
        grad_y = angle(P./circshift(P, [ 1  0  0]));
        grad_x = angle(P./circshift(P, [ 0  1  0]));
        grad_z = angle(P./circshift(P, [ 0  0  1])); clear P;
        % From which the second order derivative is calculated:
        Lap_Phase_y = grad_y - circshift(grad_y, [-1  0  0]); clear grad_y;
        Lap_Phase_x = grad_x - circshift(grad_x, [ 0 -1  0]); clear grad_x;
        Lap_Phase_z = grad_z - circshift(grad_z, [ 0  0 -1]); clear grad_x;
        % Whose summation is the Laplacian: 
        phaseLaplacianFiltered = (Lap_Phase_y+Lap_Phase_x+Lap_Phase_z)/7;
        clear Lap_Phase_y Lap_Phase_x Lap_Phase_z;
    else
        phaseLaplacianFiltered = imag(conj(P) .* ifftn(del_op .* fftn(P)));
    end
    
    %% Inverse Laplacian phase filtering:
    phaseUnwrapped = real(ifftn(del_inv.*fftn(phaseLaplacianFiltered)));
    
    %% reembedding into original matrix size:
    phaseUnwrapped_2 = zeros(Nold);
    phaseUnwrapped_2(y1, x1, z1) = phaseUnwrapped(y2, x2, z2);
    phaseUnwrapped = phaseUnwrapped_2.*msk; 
    
    phaseLaplacianFiltered_2 = zeros(Nold);
    phaseLaplacianFiltered_2(y1, x1, z1) = phaseLaplacianFiltered(y2, x2, z2);
    phaseLaplacianFiltered = phaseLaplacianFiltered_2.*msk; 
end
