rng 'default'

M = 7; % Number of observations.
N = 5; % Number of variables observed.

% Made-up data
X = rand(M,N);

% De-mean (MATLAB will de-mean inside of PCA, but I want the de-meaned values later).
X = bsxfun(@minus,X,mean(X));

% Do the PCA.
[coeff,score,latent,~,explained] = pca(X);

% Calculate eigenvalues and eigenvectors of the covariance matrix.
covarianceMatrix = cov(X);
[V,D] = eig(covarianceMatrix);

% "coeff" are the principal component vectors. These are the eigenvectors of the covariance matrix. Compare ...
coeff
V

% Multiply the original data by the principal component vectors to get the projections of the original data on the
% principal component vector space. This is also the output "score". Compare ...

dataInPrincipalComponentSpace = X*coeff
score

% The columns of X*coeff are orthogonal to each other. This is shown with ...

corrcoef(dataInPrincipalComponentSpace)

% The variances of these vectors are the eigenvalues of the covariance matrix, and are also the output "latent". Compare
% these three outputs

var(dataInPrincipalComponentSpace)'
latent
sort(diag(D),'descend')
