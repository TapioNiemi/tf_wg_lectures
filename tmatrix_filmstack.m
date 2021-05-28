%% FYS-6607 Thin-films and waveguides-course, 2021
% Tapio Niemi and Janne Simonen, 
% tapio.niemi@tuni.fi

%% Calculates reflection, r, transmission, t, reflectance, R, 
% and transmittance, T of a thin-film stack
% using the transfer matrix method.

% Formulation follows closely the approach by Yeh and Yariv
% Photonics, ISBN 13: 9780195179460, Oxford University Press
% or
% Optical Waves in Layered Media, ISBN-13: 978-0471731924
% Wiley Series in Pure and Applied Optics

%% Parameters:
% n1 = refractive index of incidence side material
% n2 = refractive index of semi-infinite substrate
% n = refractive indices of the thin-film stack [vector]
% L = thicknesses of the layers of the stack [vector] [in nanometers]
% lambda = wavelength [in nanometers]
% theta = incidence angle [in degrees]
% polarization = 'TE' or 'TM'


%% Using the function (Example 3.4):
% lambda=linspace(600,1400,100);
% for m=1:length(lambda)
% [r(m),t(m),R(m),T(m)]=tmatrix_filmstack(1,1.5,repmat([1.45 2.0],1,10),...
% repmat(1000./[1.45 2.0]/4,1,10),lambda(m),0,'TM');
% end

function [r,t,R,T]=tmatrix_filmstack(n1,n2,n,L,lambda,theta,polarization)

% Check validity of input parameters
if length(n)~=length(L)
    error('n and L must be of the same length!')
end

% refractive index vector including surrounding materials
n=[n1 n n2];

% Number of refractive indices = length of vector n
% We have N materials, N-1 interfaces and N-2 layers
N=length(n);

% Tangential  component of the incident wavevector k.
% Important: this is the same in each of the layers!
kt=2*pi/lambda*n1*sind(theta);

% Normal components of wave vectors in all layers from Pythagoras
% Notice, this is a vector of length N
kn=sqrt((2*pi/lambda*n).^2-kt^2);

% Ensure exponential decay in exit material
% This is important if total internal reflection takes place
kn=real(kn)-1j*abs(imag(kn));

% initialize the 2x2 transfer matrix
matrix=eye(2);

% This part depends on the polarization
% Reflectance and transmittance coefficients of each interface,
% vectors of length N-1

if strcmp(polarization,'TE')

    rn=(kn(1:N-1)-kn(2:N))./...
        (kn(1:N-1)+kn(2:N));
    tn=2*kn(1:N-1)./...
        (kn(1:N-1)+kn(2:N));
    
elseif strcmp(polarization,'TM')
    
    rn=(n(1:N-1).^2.*kn(2:N)-n(2:N).^2.*kn(1:N-1))./...
    (n(2:N).^2.*kn(1:N-1)+n(1:N-1).^2.*kn(2:N));
    tn=(2*n(2:N).*n(1:N-1).*kn(1:N-1))./...
        (n(2:N).^2.*kn(1:N-1)+n(1:N-1).^2.*kn(2:N));  
    
else
    error('Polarization must be TE or TM!')
end

% Generate the T and P matrices for each layer/interface
% Let's loop through all the layers, its fast enough
for M=1:N-2
    Tn=1/tn(M)*[1 rn(M);rn(M) 1];
    P=[exp(1i*kn(M+1)*L(M)) 0 ; 0 exp(-1i*kn(M+1)*L(M))];
    
    % multiply the total transfer matrix by P and T, repeat if needed
    matrix=matrix*Tn*P;
end

% We have almost the full system matrix but are missing the final
% interface T
matrix=matrix*(1/tn(N-1)*[1 rn(N-1);rn(N-1) 1]); 

% Results:
% Reflection and transmission of fields
r=(matrix(2,1)/matrix(1,1));
t=1./matrix(1,1);

% Reflectance and transmittance
R=abs(r).^2;
T=real(kn(end)./kn(1))*abs(t).^2;
end

