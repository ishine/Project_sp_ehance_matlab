doas=[-30 0 0]*pi/180; %DOA's of signals in rad.
P=[1 1 1]; %Power of incoming signals
N=2; %Number of array elements
K=1024; %Number of data snapshots
d=0.5; %Distance between elements in wavelengths
noise_var=40; %Variance of noise
r=length(doas); %Total number of signals
% Steering vector matrix. Columns will contain the steering vectors of the r signals
A=exp(-i*2*pi*d*(0:N-1)'*sin([doas(:).']));
% Signal and noise generation
sig=round(rand(r,K))*2-1; % Generate random BPSK symbols for each of the
% r signals
noise=sqrt(noise_var/2)*(randn(N,K)+i*randn(N,K)); %Uncorrelated noise
X=A*diag(sqrt(P))*sig+noise; %Generate data matrix
R=X*X'/K; %Spatial covariance matrix
%MVDR
angles = [-100:10:100];

IR=inv(R); %Inverse of covariance matrix
for k=1:length(angles)
mvdr(k)=1/(a1(:,k)'*IR*a1(:,k));
end
figure;
plot(angles,abs(mvdr)/max(abs(mvdr)),'k');hold on;
xlabel('Angle in degrees')
%Estimate DOA's using the classical beamformer
for k=1:length(angles)
Classical(k)=(a1(:,k)'*R*a1(:,k));
end
plot(angles,abs(Classical)/max(abs(Classical)),'r--');grid on;
legend('MVDR','Classical Beamformer');