%% FYS-6607 Thin-films and waveguides-course, 2021
% Tapio Niemi and Janne Simonen
% tapio.niemi@tuni.fi

%% Calculates the effective indidces of modes for a multilayer slab 
% waveguide using the transfer-matrix method. 
% Waveguide modes and their effective indices can be found when 
% M11(neff_test)=0

% n1 = refractive index of incidence side material
% n2 = refractive index after the thin film stack
% n = refractive indices of the layers of the stack, vector
% L = thicknesses of the layers of the stack, vector [in nanometers]
% neff_test = (test) effective index of the mode. Scan this to find M11=0.
% polarization = 'TE' or 'TM'


%% Using the function: 
% Reference results (Yariv, Photonics)
% Refractive indices: cladding1=1, core=2, cladding2=1.7
% thickness=wavelength=800 nm

% TM : fminsearch(@(neff)
% abs(M11_multilayer_waveguide(1,1.7,2,800e-9,800e-9,neff,'TM')),1.96) 
% =1.9513

% OR

% fzero(@(neff)
% real(M11_multilayer_waveguide(1,1.7,2,800e-9,800e-9,neff,'TM')),1.96)
% =1.9513

% TE : fminsearch(@(neff)
% abs(M11_multilayer_waveguide(1,1.7,2,800e-9,800e-9,neff,'TE')),1.96) 
% =1.9594

function M11=M11_multilayer_waveguide(n1,n2,n,L,lambda,neff_test,polarization)

% Check validity of input parameters
if length(n)~=length(L)
    error('n and L must be of the same length!')
end

% Refractive index vector including surrounding materials
n=[n1 n n2];

% Number of refractive indices = length of vector n
% (so we have N materials, N-1 interfaces and N-2 layers)
N=length(n);

% Tangential (y) component of the incident wavevector k.
% Important: this is the same in each of the layers!
kt=2*pi./lambda.*neff_test;

% Normal (x) components of wave vectors in all layers from Pythagoras
% Notice, this is a vector of length N
kn=sqrt((2*pi/lambda*n).^2-kt.^2);

% Ensure exponential decay in claddings
% In this case imaginary part should be negative
kn(1)=real(kn(1))-1j*abs(imag(kn(1)));
kn(end)=real(kn(end))-1j*abs(imag(kn(end)));

% Initialize the 2x2 transfer matrix with ones
matrix=eye(2);

% This part depends on polarization.
% Reflection and transmission coefficients of each interface
if strcmp(polarization,'TE')
    
    r=(kn(1:N-1)-kn(2:N))./...
        (kn(1:N-1)+kn(2:N));
    t=2*kn(1:N-1)./...
        (kn(1:N-1)+kn(2:N));
    
    elseif strcmp(polarization,'TM')
    
    r=(n(1:N-1).^2.*kn(2:N)-n(2:N).^2.*kn(1:N-1))./...
        (n(2:N).^2.*kn(1:N-1)+n(1:N-1).^2.*kn(2:N));
    t=(2*n(2:N).*n(1:N-1).*kn(1:N-1))./...
        (n(2:N).^2.*kn(1:N-1)+n(1:N-1).^2.*kn(2:N));  
    
else
    
    error('Polarization must be TE or TM!')
    
end

% Generate the T and P matrices for each layer/interface
% Lets loop through all the layers, its fast enough
for M=1:N-2
    T=1/t(M)*[1 r(M);r(M) 1];
    P=[exp(1i*kn(M+1)*L(M)) 0 ; 0 exp(-1i*kn(M+1)*L(M))];
    
    % Update the total transfer matrix 
    matrix=matrix*T*P;
    
end

% We have almost the full system matrix 
% but are missing the final interface T
matrix=matrix*(1/t(N-1)*[1 r(N-1);r(N-1) 1]); 
M11=matrix(1,1);
end
