%% FYS-6607 Thin-films and waveguides-course, 2021
% Tapio Niemi 
% tapio.niemi@tuni.fi

%% Calculates band structure of periodic layer stack
% using transfer matrices

% Formulation follows closely the approach by Yeh and Yariv
% Photonics, ISBN 13: 9780195179460, Oxford University Press
% or
% Optical Waves in Layered Media, ISBN-13: 978-0471731924
% Wiley Series in Pure and Applied Optics

%% Parameters:
% n = refractive indices of the thin-film stack [vector]
% L = thicknesses of the layers of the stack [vector] [in nanometers]
% lambda = wavelength [in nanometers]
% theta = incidence angle [in degrees]
% polarization = 'TE' or 'TM'

function K=tmatrix_filmstack_periodic(n,L,lambda,theta,polarization)

% Check validity of input parameters
if length(n)~=length(L)
    error('n and L must be of the same length!')
end

% refractive index vector including surrounding materials
% this will effectivelly make the structure periodic
n=[n(end) n n(1)];

% Number of refractive indices = length of vector n
% (so we have N materials, N-1 interfaces and N-2 layers)
N=length(n);

% Tangential  component of the incident wavevector k.
% Important: this is the same in each of the layers!
% Remember that sin() wants radians but theta is in degrees
kt=2*pi/lambda*n(2)*sin(theta*pi/180);

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
    r=(kn(1:N-1)-kn(2:N))./(kn(1:N-1)+kn(2:N));
    t=2*kn(1:N-1)./(kn(1:N-1)+kn(2:N));
    
elseif strcmp(polarization,'TM')
    r=(n(1:N-1).^2.*kn(2:N)-n(2:N).^2.*kn(1:N-1))./(n(2:N).^2.*kn(1:N-1)+n(1:N-1).^2.*kn(2:N));
    t=(2*n(2:N).*n(1:N-1).*kn(1:N-1))./(n(2:N).^2.*kn(1:N-1)+n(1:N-1).^2.*kn(2:N));  
else
    error('Polarization must be TE or TM!')
end

% Generate the T and P matrices for each layer/interface
% Let's loop through all the layers, its fast enough

for M=1:N-2
    T=1/t(M)*[1 r(M);r(M) 1];
    P=[exp(1i*kn(M+1)*L(M)) 0 ; 0 exp(-1i*kn(M+1)*L(M))];
    % multiply the total transfer matrix by P and T, repeat if needed
    matrix=matrix*T*P;
end

% K=1/L*acos(0.5*(A+D))
K=1/(sum(L))*acos(0.5*(matrix(1,1)+matrix(2,2)));

end

