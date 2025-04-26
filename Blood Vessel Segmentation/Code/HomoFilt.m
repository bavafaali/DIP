function [Hf] = HomoFilt(sigma,I)
    Ihf = im2double(I);
    Ihf  = log(1 + Ihf );
    M = 2*size(Ihf ,1) + 1;
    N = 2*size(Ihf ,2) + 1;
    [X, Y] = meshgrid(1:N,1:M);
    centerX = ceil(N/2); 
    centerY = ceil(M/2); 
    gaussianNumerator = (X - centerX).^2 + (Y - centerY).^2;
    H = exp(-gaussianNumerator./(2*sigma.^2));
    H = 1 - H; 
    H = fftshift(H);
    If = fft2(Ihf , M, N);
    Iout = real(ifft2(H.*If));
    Iout = Iout(1:size(Ihf ,1),1:size(Ihf ,2));
    Ihmf = im2uint8(exp(Iout) - 1);
    Hf = Ihmf;
end