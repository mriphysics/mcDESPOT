%% Load SPGR data.

SPGRFA2_Image = csmtData{1};
SPGRFA4_Image = csmtData{2};
SPGRFA6_Image = csmtData{3};
SPGRFA8_Image = csmtData{4};
SPGRFA10_Image = csmtData{5};
SPGRFA12_Image = csmtData{6};
SPGRFA14_Image = csmtData{7};
SPGRFA16_Image = csmtData{8};
SPGRFA18_Image = csmtData{9};
SPGRFA20_Image = csmtData{10};

%% Visualise data.

figure(1)
Signal = zeros(1,10);
for ii = 1:10
    Slice = csmtData{ii};
    subplot(3,5,ii); imagesc(squeeze(abs(Slice(:,:,60))));
    axis square; get(gca,'XTickLabels'); set(gca,'XTickLabels',[]);
    get(gca,'YTickLabels'); set(gca,'YTickLabels',[]);
    Signal(ii) = squeeze(abs(Slice(140,73,81)));
end

figure(1); subplot(3,5,[11 12 13 14 15])
plot(2:2:20,Signal)
xlabel('FA (deg)'); ylabel('Signal')

%% Calculate SNR.

alldata = cell2mat(csmtData);
alldata = reshape(alldata,[304 145 10 127]);
alldata = permute(alldata,[1 2 4 3]);

mask = (sum(abs(alldata),4)>5)&(sum(abs(alldata),4)<10);
mask_rep = repmat(mask,[1 1 1 10]);
% Drop some volumes?
% mask_rep(:,:,:,1)=0;

noisedist = alldata(mask_rep);


figure(2);
clf

subplot(131)
histogram(abs(noisedist),linspace(0,5,100))
subplot(132)
histogram(real(noisedist),linspace(-5,5,100))
hold on
histogram(imag(noisedist),linspace(-5,5,100))
subplot(133)
imagesc(squeeze(abs(mask(:,:,60))))

fprintf(1,'Noise using SNR_mean = %1.3f, SNR_stdv = %1.3f, std(real)=%1.3f, std(imag)=%1.3f\n',...
    sqrt(2/pi)*mean(abs(noisedist)),sqrt(2/(4-pi))*std(abs(noisedist)),std(real(noisedist)),std(imag(noisedist)));

% Crude signal estimate.
mask = (sum(abs(alldata),4)>400);
mask(175:end,:,:)=0;
mask_rep = repmat(mask,[1 1 1 10]);
% mask_rep(:,:,:,1)=0;
signaldist = alldata(mask_rep);

% Crude signal estimate v2.
% imagesc(squeeze(abs(SPGRFA14_Image(140,:,:))))
% mask = roipoly;

figure(3);
clf

subplot(121)
histogram(abs(signaldist))

subplot(122)
imagesc(squeeze(abs(mask(:,:,60))))
fprintf(1,'signal = %1.3f\n',mean(abs(signaldist)))

%% Initial test - same ROI.

imagesc(squeeze(abs(SPGRFA14_Image(140,:,:))));
roi = roipoly;

fa = [2:2:20];
snr=[];
for ii=1:10
    eval(sprintf('tmp=squeeze(SPGRFA%d_Image(140,:,:));',fa(ii)));
    tmp2=abs(tmp(roi));
    snr(ii)= mean((tmp2))/std((tmp2));
end