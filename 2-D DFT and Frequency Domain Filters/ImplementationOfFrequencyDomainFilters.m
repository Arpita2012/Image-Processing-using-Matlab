

% Calculation of 2-D Discrete Fourier Transform 
% [Conversion from spatial domain to frequency domain]
% And implementation of frequency domain filters
clear all
%read input image
Image=imread('cameraman.tif');
%Image=imread('ex1.jpg'); //Given Image 
Img=double(Image);
[r,c]=size(Img);

r1=2*r;
c1=2*c;

pim=zeros((r1),(c1));
paddedImg=zeros((r1),(c1));

%padding
for i=1:r
    for j=1:c
   pim(i,j)=Img(i,j);
    end
end

%Center the transform
%[ multiply the input image by (-1)x+y to center the 
%transform to u = M/2 and v = N/2 (if M and N are 
%even numbers, then the shifted coordinates will be integers)]
for i=1:r
    for j=1:c
   paddedImg(i,j)=pim(i,j)*((-1)^(i+j));
    end
end


%2D fft
fourierTransformedImage=fft2(paddedImg);

n=1; %order for butterworth filter
thresh=30; % cutoff radius in frequency domain for filters





%% GAUSSIAN LOW PASS FILTER
 gaussianLowPassFilteredImage=glp(fourierTransformedImage,thresh); % gaussian low pass filter
% him=blpf(fim,thresh,n); % butterworth low pass filter

% % function calls for high pass filters
% him=ghp(fim,thresh); % gaussian low pass filter
%  him=bhp(fim,thresh,n);  %butterworth high pass filter


%inverse 2D fft
 inverseFourierTransform_GaussianLowPassFilteredImage=ifft2(gaussianLowPassFilteredImage);
 
for i=1:r1
    for j=1:c1
   inverseFourierTransform_GaussianLowPassFilteredImage(i,j)=inverseFourierTransform_GaussianLowPassFilteredImage(i,j)*((-1)^(i+j));
    end
end


% removing the padding
for i=1:r
    for j=1:c
   inverseFourierTransform_GaussianLowPassFilteredImageZeroPadding(i,j)=inverseFourierTransform_GaussianLowPassFilteredImage(i,j);
    end
end

% retaining the real parts of the matrix
inverseFourierTransform_GaussianLowPassFilteredImageZeroPadding=real(inverseFourierTransform_GaussianLowPassFilteredImageZeroPadding);
inverseFourierTransform_GaussianLowPassFilteredImageZeroPadding=uint8(inverseFourierTransform_GaussianLowPassFilteredImageZeroPadding);

%figure, imshow(inverseFourierTransform_GaussianLowPassFilteredImageZeroPadding);title('Using Gaussian Low Pass Filter');

figure;


 subplot(2,3,1);imshow(Image);title('Original image');
 subplot(2,3,2);imshow(uint8(paddedImg));title('1. Padding');
 subplot(2,3,3);imshow(uint8(fourierTransformedImage));title('2.Compute F(u, v), the 2-D DFT of the image (1)');
 subplot(2,3,4);imshow(uint8(gaussianLowPassFilteredImage));title('3. Multiply F(u, v) by a filter function H(u, v)');
 subplot(2,3,5);imshow(uint8(inverseFourierTransform_GaussianLowPassFilteredImage));title('4. Compute the inverse 2-D DFT of the result in (3)');
  text(size(inverseFourierTransform_GaussianLowPassFilteredImage,2),size(inverseFourierTransform_GaussianLowPassFilteredImage,1)+15, ...
    '5. Displaying Only real Part of result in (4) ', ...
    'FontSize',10,'HorizontalAlignment','Right');
 subplot(2,3,6);imshow(uint8(inverseFourierTransform_GaussianLowPassFilteredImageZeroPadding));title('6.Resultant Image [Padding Removed] ');
 text(size(inverseFourierTransform_GaussianLowPassFilteredImageZeroPadding,2),size(inverseFourierTransform_GaussianLowPassFilteredImageZeroPadding,1)+15, ...
    'Using Gaussian Low Pass Filter', ...
    'FontSize',10,'HorizontalAlignment','Right');







%% BUTTERWORTH LOW PASS FILTER
butterworthLowPassFilteredImage=blpf(fourierTransformedImage,thresh,n); % butterworth low pass filter

% % function calls for high pass filters
% him=ghp(fim,thresh); % gaussian low pass filter
%  him=bhp(fim,thresh,n);  %butterworth high pass filter





%inverse 2D fft
 inverseFourierTransform_ButterworthLowPassFilteredImage=ifft2(butterworthLowPassFilteredImage);
 
for i=1:r1
    for j=1:c1
   inverseFourierTransform_ButterworthLowPassFilteredImage(i,j)=inverseFourierTransform_ButterworthLowPassFilteredImage(i,j)*((-1)^(i+j));
    end
end


% removing the padding
for i=1:r
    for j=1:c
   inverseFourierTransform_ButterworthLowPassFilteredImageZeroPadding(i,j)=inverseFourierTransform_ButterworthLowPassFilteredImage(i,j);
    end
end

% retaining the real parts of the matrix
inverseFourierTransform_ButterworthLowPassFilteredImageZeroPadding=real(inverseFourierTransform_ButterworthLowPassFilteredImageZeroPadding);
inverseFourierTransform_ButterworthLowPassFilteredImageZeroPadding=uint8(inverseFourierTransform_ButterworthLowPassFilteredImageZeroPadding);

%figure, imshow(inverseFourierTransform_ButterworthLowPassFilteredImageZeroPadding);title('Using Butterworth Low Pass Filter');

figure;


 subplot(2,3,1);imshow(Image);title('Original image');
 subplot(2,3,2);imshow(uint8(paddedImg));title('1. Padding');
 subplot(2,3,3);imshow(uint8(fourierTransformedImage));title('2.Compute F(u, v), the 2-D DFT of the image (1)');
 subplot(2,3,4);imshow(uint8(butterworthLowPassFilteredImage));title('3. Multiply F(u, v) by a filter function H(u, v)');
 subplot(2,3,5);imshow(uint8(inverseFourierTransform_ButterworthLowPassFilteredImage));title('4. Compute the inverse 2-D DFT of the result in (3)');
  text(size(inverseFourierTransform_ButterworthLowPassFilteredImage,2),size(inverseFourierTransform_ButterworthLowPassFilteredImage,1)+15, ...
    '5. Displaying Only real Part of result in (4) ', ...
    'FontSize',10,'HorizontalAlignment','Right');
 subplot(2,3,6);imshow(uint8(inverseFourierTransform_ButterworthLowPassFilteredImageZeroPadding));title('6.Resultant Image [Padding Removed] ');
 text(size(inverseFourierTransform_ButterworthLowPassFilteredImageZeroPadding,2),size(inverseFourierTransform_ButterworthLowPassFilteredImageZeroPadding,1)+15, ...
    'Using Butterworth Low Pass Filter', ...
    'FontSize',10,'HorizontalAlignment','Right');




%% GAUSSIAN HIGH PASS FILTER
 gaussianHighPassFilteredImage=ghp(fourierTransformedImage,thresh); % gaussian High pass filter

 
%inverse 2D fft
 inverseFourierTransform_GaussianHighPassFilteredImage=ifft2(gaussianHighPassFilteredImage);
 
for i=1:r1
    for j=1:c1
   inverseFourierTransform_GaussianHighPassFilteredImage(i,j)=inverseFourierTransform_GaussianHighPassFilteredImage(i,j)*((-1)^(i+j));
    end
end
 
 
% removing the padding
for i=1:r
    for j=1:c
   inverseFourierTransform_GaussianHighPassFilteredImageZeroPadding(i,j)=inverseFourierTransform_GaussianHighPassFilteredImage(i,j);
    end
end
 
% retaining the real parts of the matrix
inverseFourierTransform_GaussianHighPassFilteredImageZeroPadding=real(inverseFourierTransform_GaussianHighPassFilteredImageZeroPadding);
inverseFourierTransform_GaussianHighPassFilteredImageZeroPadding=uint8(inverseFourierTransform_GaussianHighPassFilteredImageZeroPadding);
 
%figure, imshow(inverseFourierTransform_GaussianHighPassFilteredImageZeroPadding);title('Using Gaussian High Pass Filter');
 
figure;
 
 
 subplot(2,3,1);imshow(Image);title('Original image');
 subplot(2,3,2);imshow(uint8(paddedImg));title('1. Padding');
 subplot(2,3,3);imshow(uint8(fourierTransformedImage));title('2.Compute F(u, v), the 2-D DFT of the image (1)');
 subplot(2,3,4);imshow(uint8(gaussianHighPassFilteredImage));title('3. Multiply F(u, v) by a filter function H(u, v)');
 subplot(2,3,5);imshow(uint8(inverseFourierTransform_GaussianHighPassFilteredImage));title('4. Compute the inverse 2-D DFT of the result in (3)');
  text(size(inverseFourierTransform_GaussianHighPassFilteredImage,2),size(inverseFourierTransform_GaussianHighPassFilteredImage,1)+15, ...
    '5. Displaying Only real Part of result in (4) ', ...
    'FontSize',10,'HorizontalAlignment','Right');
 subplot(2,3,6);imshow(uint8(inverseFourierTransform_GaussianHighPassFilteredImageZeroPadding));title('6.Resultant Image [Padding Removed] ');
 text(size(inverseFourierTransform_GaussianHighPassFilteredImageZeroPadding,2),size(inverseFourierTransform_GaussianHighPassFilteredImageZeroPadding,1)+15, ...
    'Using Gaussian High Pass Filter', ...
    'FontSize',10,'HorizontalAlignment','Right');
 



%% BUTTERWORTH HIGH PASS FILTER
butterworthHighPassFilteredImage=bhp(fourierTransformedImage,thresh,n); % butterworth High pass filter
 
   
%inverse 2D fft
 inverseFourierTransform_ButterworthHighPassFilteredImage=ifft2(butterworthHighPassFilteredImage);
 
for i=1:r1
    for j=1:c1
   inverseFourierTransform_ButterworthHighPassFilteredImage(i,j)=inverseFourierTransform_ButterworthHighPassFilteredImage(i,j)*((-1)^(i+j));
    end
end
 
 
% removing the padding
for i=1:r
    for j=1:c
   inverseFourierTransform_ButterworthHighPassFilteredImageZeroPadding(i,j)=inverseFourierTransform_ButterworthHighPassFilteredImage(i,j);
    end
end
 
% retaining the real parts of the matrix
inverseFourierTransform_ButterworthHighPassFilteredImageZeroPadding=real(inverseFourierTransform_ButterworthHighPassFilteredImageZeroPadding);
inverseFourierTransform_ButterworthHighPassFilteredImageZeroPadding=uint8(inverseFourierTransform_ButterworthHighPassFilteredImageZeroPadding);
 
%figure, imshow(inverseFourierTransform_ButterworthHighPassFilteredImageZeroPadding);title('Using Butterworth High Pass Filter');
 
figure;
 
 
 subplot(2,3,1);imshow(Image);title('Original image');
 subplot(2,3,2);imshow(uint8(paddedImg));title('1. Padding');
 subplot(2,3,3);imshow(uint8(fourierTransformedImage));title('2.Compute F(u, v), the 2-D DFT of the image (1)');
 subplot(2,3,4);imshow(uint8(butterworthHighPassFilteredImage));title('3. Multiply F(u, v) by a filter function H(u, v)');
 subplot(2,3,5);imshow(uint8(inverseFourierTransform_ButterworthHighPassFilteredImage));title('4. Compute the inverse 2-D DFT of the result in (3)');
  text(size(inverseFourierTransform_ButterworthHighPassFilteredImage,2),size(inverseFourierTransform_ButterworthHighPassFilteredImage,1)+15, ...
    '5. Displaying Only real Part of result in (4) ', ...
    'FontSize',10,'HorizontalAlignment','Right');
 subplot(2,3,6);imshow(uint8(inverseFourierTransform_ButterworthHighPassFilteredImageZeroPadding));title('6.Resultant Image [Padding Removed] ');
 text(size(inverseFourierTransform_ButterworthHighPassFilteredImageZeroPadding,2),size(inverseFourierTransform_ButterworthHighPassFilteredImageZeroPadding,1)+15, ...
    'Using Butterworth High Pass Filter', ...
    'FontSize',10,'HorizontalAlignment','Right');
 


 


%%  IDEAL LOW PASS FILTER
 idealLowPassFilteredImage=ilp(fourierTransformedImage,thresh); % ideal low pass filter
% him=blpf(fim,thresh,n); % butterworth low pass filter
 
% % function calls for high pass filters
% him=ghp(fim,thresh); % ideal low pass filter
%  him=bhp(fim,thresh,n);  %butterworth high pass filter
 
 
%inverse 2D fft
 inverseFourierTransform_IdealLowPassFilteredImage=ifft2(idealLowPassFilteredImage);
 
for i=1:r1
    for j=1:c1
   inverseFourierTransform_IdealLowPassFilteredImage(i,j)=inverseFourierTransform_IdealLowPassFilteredImage(i,j)*((-1)^(i+j));
    end
end
 
 
% removing the padding
for i=1:r
    for j=1:c
   inverseFourierTransform_IdealLowPassFilteredImageZeroPadding(i,j)=inverseFourierTransform_IdealLowPassFilteredImage(i,j);
    end
end
 
% retaining the real parts of the matrix
inverseFourierTransform_IdealLowPassFilteredImageZeroPadding=real(inverseFourierTransform_IdealLowPassFilteredImageZeroPadding);
inverseFourierTransform_IdealLowPassFilteredImageZeroPadding=uint8(inverseFourierTransform_IdealLowPassFilteredImageZeroPadding);
 
%figure, imshow(inverseFourierTransform_IdealLowPassFilteredImageZeroPadding);title('Using Ideal Low Pass Filter');
 
figure;
 
 
 subplot(2,3,1);imshow(Image);title('Original image');
 subplot(2,3,2);imshow(uint8(paddedImg));title('1. Padding');
 subplot(2,3,3);imshow(uint8(fourierTransformedImage));title('2.Compute F(u, v), the 2-D DFT of the image (1)');
 subplot(2,3,4);imshow(uint8(idealLowPassFilteredImage));title('3. Multiply F(u, v) by a filter function H(u, v)');
 subplot(2,3,5);imshow(uint8(inverseFourierTransform_IdealLowPassFilteredImage));title('4. Compute the inverse 2-D DFT of the result in (3)');
  text(size(inverseFourierTransform_IdealLowPassFilteredImage,2),size(inverseFourierTransform_IdealLowPassFilteredImage,1)+15, ...
    '5. Displaying Only real Part of result in (4) ', ...
    'FontSize',10,'HorizontalAlignment','Right');
 subplot(2,3,6);imshow(uint8(inverseFourierTransform_IdealLowPassFilteredImageZeroPadding));title('6.Resultant Image [Padding Removed] ');
 text(size(inverseFourierTransform_IdealLowPassFilteredImageZeroPadding,2),size(inverseFourierTransform_IdealLowPassFilteredImageZeroPadding,1)+15, ...
    'Using Ideal Low Pass Filter', ...
    'FontSize',10,'HorizontalAlignment','Right');
 
 


%% IDEAL HIGH PASS FILTER
 idealHighPassFilteredImage=ihp(fourierTransformedImage,thresh); % ideal high pass filter
% him=blpf(fim,thresh,n); % butterworth high pass filter
 
% % function calls for high pass filters
% him=ghp(fim,thresh); % ideal high pass filter
%  him=bhp(fim,thresh,n);  %butterworth high pass filter
 
 
%inverse 2D fft
 inverseFourierTransform_IdealHighPassFilteredImage=ifft2(idealHighPassFilteredImage);
 
for i=1:r1
    for j=1:c1
   inverseFourierTransform_IdealHighPassFilteredImage(i,j)=inverseFourierTransform_IdealHighPassFilteredImage(i,j)*((-1)^(i+j));
    end
end
 
 
% removing the padding
for i=1:r
    for j=1:c
   inverseFourierTransform_IdealHighPassFilteredImageZeroPadding(i,j)=inverseFourierTransform_IdealHighPassFilteredImage(i,j);
    end
end
 
% retaining the ral parts of the matrix
inverseFourierTransform_IdealHighPassFilteredImageZeroPadding=real(inverseFourierTransform_IdealHighPassFilteredImageZeroPadding);
inverseFourierTransform_IdealHighPassFilteredImageZeroPadding=uint8(inverseFourierTransform_IdealHighPassFilteredImageZeroPadding);
 
%figure, imshow(inverseFourierTransform_IdealHighPassFilteredImageZeroPadding);title('Using Ideal High Pass Filter');
 
figure;
 
 
 subplot(2,3,1);imshow(Image);title('Original image');
 subplot(2,3,2);imshow(uint8(paddedImg));title('1. Padding');
 subplot(2,3,3);imshow(uint8(fourierTransformedImage));title('2.Compute F(u, v), the 2-D DFT of the image (1)');
 subplot(2,3,4);imshow(uint8(idealHighPassFilteredImage));title('3. Multiply F(u, v) by a filter function H(u, v)');
 subplot(2,3,5);imshow(uint8(inverseFourierTransform_IdealHighPassFilteredImage));title('4. Compute the inverse 2-D DFT of the result in (3)');
  text(size(inverseFourierTransform_IdealHighPassFilteredImage,2),size(inverseFourierTransform_IdealHighPassFilteredImage,1)+15, ...
    '5. Displaying Only real Part of result in (4) ', ...
    'FontSize',10,'HorizontalAlignment','Right');
 subplot(2,3,6);imshow(uint8(inverseFourierTransform_IdealHighPassFilteredImageZeroPadding));title('6.Resultant Image [Padding Removed] ');
 text(size(inverseFourierTransform_IdealHighPassFilteredImageZeroPadding,2),size(inverseFourierTransform_IdealHighPassFilteredImageZeroPadding,1)+15, ...
    'Using Ideal High Pass Filter', ...
    'FontSize',10,'HorizontalAlignment','Right');
 
 


%% SUMMARY


figure;
 
 
 subplot(2,3,1);imshow(inverseFourierTransform_IdealLowPassFilteredImageZeroPadding);title('Using Ideal Low Pass Filter');
 subplot(2,3,2);imshow(inverseFourierTransform_GaussianLowPassFilteredImageZeroPadding);title('Using Gaussian Low Pass Filter');
 subplot(2,3,3);imshow(inverseFourierTransform_ButterworthLowPassFilteredImageZeroPadding);title('Using Butterworth Low Pass Filter');
 subplot(2,3,4);imshow(inverseFourierTransform_IdealHighPassFilteredImageZeroPadding);title('Using Ideal High Pass Filter');
 subplot(2,3,5);imshow(inverseFourierTransform_GaussianHighPassFilteredImageZeroPadding);title('Using Gaussian High Pass Filter');
 subplot(2,3,6);imshow(inverseFourierTransform_ButterworthHighPassFilteredImageZeroPadding);title('Using Butterworth High Pass Filter');
 


