%% FYS-6607 Thin-films and waveguides-course, 2021
% Tapio Niemi  
% tapio.niemi@tuni.fi

%% Numerical evaluation of Fraunhofer and Fresnel diffractions in 1D
% y0 = axis for input electric field [m]
% E_field = complex electric field
% z = diffraction distance [m]
% lambda = wavelength [m]
% model = 'Fresnel' or ' Fraunhofer'

function [y E_Out] = diffraction_FFT(y0,E_field,z,lambda,model)

% Scaling of the diffraction axis
dydot=abs(y0(2)-y0(1));     % Sampling distance of the field
N=length(y0);                % Number of points
k=2*pi/lambda;

% New y-axis for the transformed field
y = lambda*z/(dydot*N)*linspace(-N/2,N/2-1,N);
Complex_amplitude=sqrt(j/z/lambda)*exp(-j*k*z)*exp(-j*k/(2*z)*y.^2);

if strcmp(model,'Fresnel')

    E_Out=dydot*Complex_amplitude.*...
        fftshift(fft(ifftshift(exp(-j*k/(2*z)*y0.^2).*E_field)));

elseif strcmp(model,'Fraunhofer')

    E_Out=dydot*Complex_amplitude.*fftshift(fft(ifftshift(E_field)));

else
    error('must be Fresnel or Fraunhofer!')
end

% FFT- fftshift moves and flips the result around center
% This should be used to keep the phase of the field in order
% Undocumented "feature" of MATLAB
% For normal spectral analysis fftshift(fft((E_field)) works just fine

end

