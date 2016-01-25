function [] = testForensic2D(imageScript)
% Take image
x = double(rgb2gray(imread(imageScript)));
x = x - min(x(:));
x = x/max(x(:)); 

%Choose a random alpha
alpha = randn(3);
alpha(2,2) = 0;

% Choose N and sigma 
sigma = 0.0075; %Standard Deviation

% Set p0 to the range of the reciprocal of the signal y
p0 = 1/(max(x(:))+0.01);

%Set Y as in eq.14

%Set h to be a binomial low pass filter of size (Nh x Nh)
h = [1 2 1; 2 4 2; 1 2 1]./16;

n = 0;
[Y,Z] = size(x);
while(1),
    % Expectation step
    R = zeros(size(x)); p = R; w = R; E = R;
    for i = 2:Y-1,
        for j = 2:Z-1,
            R(i,j) = abs(x(i,j) - sum(sum(alpha.*x(i-1:i+1,j-1:j+1))));
        end
    end
    for i = 2:Y-1,
        for j = 2:Z-1,
            p(i,j) = (1/(sigma*sqrt(2*pi))) * exp((-R(i,j)^2) /  (2 * sigma^2)); %Conditional Probability
            w(i,j) = p(i,j)/(p(i,j) + p0); %Posterior Probability
            E(i,j) = w(i,j)*(R(i,j)^2);
        end
    end
    %Maximization step
    sigma = sqrt(sum(E(:))/sum(w(:))); % new variance
    rightMat = [];B = zeros(9,9);col = 1;
    for s = 1:3,
        for t = 1:3,
            leftMat = [];
            for u = 1:3,
                for v = 1:3,
                    sum_lhs = 0; sum_rhs = 0;
                    for i = 2:Y-2,
                        for j = 2:Z-2,
                            sum_lhs = sum_lhs + w(i,j)*x(i+s-2,j+t-2)*x(i+u-2,j+v-2);
                            sum_rhs = sum_rhs + w(i,j)*x(i+s-2,j+t-2)*x(i,j);
                        end
                    end
                    leftMat = [leftMat; sum_lhs];
                end
            end
            rightMat = [rightMat sum_rhs];
            B(col,:) = leftMat;
            col = col+1;
        end
    end
    B_inv = (B)^(-1);
    alphaUV = rightMat * B_inv;
    alpha_new = [alphaUV(1:3);alphaUV(4:6);alphaUV(7:9)]; 
    alpha_new(2,2) = 0;
    n = n + 1; fprintf(strcat(num2str(n),' '));
    if norm(alpha-alpha_new,2) < 0.01,
        break
    else
        alpha = alpha_new;
    end
end

H = padarray(2,[2 2]) - fspecial('gaussian' ,[5 5],2);
P = fftshift(abs(fft2(p)));
Ph = imfilter(P,H);
Pg = ((Ph./max(abs(Ph(:)))).^4) * max(abs(Ph(:)));

[J,rectTampered] = imcrop(x);
[J,rectTestOrig] = imcrop(x);

I2 = imcrop(p,floor(rectTampered));
I3 = imcrop(p,floor(rectTestOrig));
P = fftshift(abs(fft2(I2)));
Ph = imfilter(P,H);
Pg = ((Ph./max(abs(Ph(:)))).^4) * max(abs(Ph(:)));
figure;imshow(Pg,[]);
P = fftshift(abs(fft2(I3)));
Ph = imfilter(P,H);
Pg = ((Ph./max(abs(Ph(:)))).^4) * max(abs(Ph(:)));
figure;imshow(Pg,[]);

end


% Results 1
% I2 = imcrop(p,[840,180,160,100]);
% I3 = imcrop(p,[1000,180,160,100]);
% P = fftshift(abs(fft2(I2)));
% Ph = imfilter(P,H);
% figure;imshow(Ph,[]);
% Pg = ((Ph./max(abs(Ph(:)))).^4) * max(abs(Ph(:)));
% figure;imshow(Pg,[]);
% P = fftshift(abs(fft2(I3)));
% Ph = imfilter(P,H);
% figure;imshow(Ph,[]);
% Pg = ((Ph./max(abs(Ph(:)))).^4) * max(abs(Ph(:)));
% figure;imshow(Pg,[]);

% alpha = (Results 1)
% 
%    1.0e-12 *
% 
%    -0.2558    0.3411         0
%     0.2274         0   -0.1634
%    -0.1421    0.1705    0.0568


% Results 2
% I2 = imcrop(p,[1,840,470,260]);
% I3 = imcrop(p,[2900,1310,470,260]);
% P = fftshift(abs(fft2(I2)));
% Ph = imfilter(P,H);
% figure;imshow(Ph,[]);
% Pg = ((Ph./max(abs(Ph(:)))).^4) * max(abs(Ph(:)));
% figure;imshow(Pg,[]);
% P = fftshift(abs(fft2(I3)));
% Ph = imfilter(P,H);
% figure;imshow(Ph,[]);
% Pg = ((Ph./max(abs(Ph(:)))).^4) * max(abs(Ph(:)));
% figure;imshow(Pg,[]);

% alpha = (for results 2)
% 
%    1.0e-12 *
% 
%    -0.3411         0   -0.6821
%    -0.5684         0    0.2842
%    -0.0284    0.4547   -0.4547

% Results 3
% I2 = imcrop(p,[1120,1140,1440,692]);
% I3 = imcrop(p,[400,1400,1440,692]);
% P = fftshift(abs(fft2(I2)));
% Ph = imfilter(P,H);
% figure;imshow(Ph,[]);
% Pg = ((Ph./max(abs(Ph(:)))).^4) * max(abs(Ph(:)));
% figure;imshow(Pg,[]);
% P = fftshift(abs(fft2(I3)));
% Ph = imfilter(P,H);
% figure;imshow(Ph,[]);
% Pg = ((Ph./max(abs(Ph(:)))).^4) * max(abs(Ph(:)));
% figure;imshow(Pg,[]);

% alpha = (result 3)
% 
%    1.0e-12 *
% 
%    -0.1990    0.3411   -0.0284
%     0.0568         0   -0.1421
%     0.0107    0.0853   -0.0853