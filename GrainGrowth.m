function GrainGrowth(Nx, Ny, dx, dy, end_time, time_step)
%   Phase-field simulation code for grain growth in 2D, using the Allen-Cahn
%   equation (for non-conserved order parameters).
%   
%   The model is developed by Fan & Chen (http://www.mmm.psu.edu/DNFan1997Actamater_Graingrowth1phase2D.pdf)
%   
%   Function call: GrainGrowth(Nx, Ny, dx, dy, end_time, time_step)
%
%   The current simulation uses 10 order parameters for 10 different grain
%   orientations.
%   
%   Output: PNG snapshot of the microstructures at 'time_step' intervals
%   Nx : Grid dimension in x direction
%   Ny : Grid dimension in y direction
%   dx : Grid spacing in x direction
%   dy : Grid spacing in y direction
%   end_time : End time of the simulation run
%   time_step : the time interval for output of simulation data
%
%   Copyright (C) 2021  Abhinav Roy
%
%   This program is free software: you can redistribute it and/or modify
%   it under the terms of the GNU General Public License as published by
%   the Free Software Foundation, either version 3 of the License, or
%   any later version.
%   This program is distributed in the hope that it will be useful,
%   but WITHOUT ANY WARRANTY; without even the implied warranty of
%   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
%   GNU General Public License for more details.
%   You should have received a copy of the GNU General Public License
%   along with this program. If not, see <https://www.gnu.org/licenses/>.

tic
if nargin < 6
    error("Please enter the required number of arguments. For more on in function, use help GrainGrowth2D");
end
disp('The code execution has commenced');
axis("square");
more off;
%--------------------------------------------------------------------------------------------------
% SIMULATION PARAMETERS
halfNx = Nx/2;
halfNy = Ny/2;
start_time = 1;  
%--------------------------------------------------------------------------------------------------
% Defining the initial profile of the order parameters
phi = unidrnd(10,Nx,Ny);
% Defining the ten order parameters
eta1 = zeros(Nx,Ny);
eta2 = zeros(Nx,Ny);
eta3 = zeros(Nx,Ny);
eta4 = zeros(Nx,Ny);
eta5 = zeros(Nx,Ny);
eta6 = zeros(Nx,Ny);
eta7 = zeros(Nx,Ny);
eta8 = zeros(Nx,Ny);
eta9 = zeros(Nx,Ny);
eta10 = zeros(Nx,Ny);
%--------------------------------------------------------------------------------------------------
% Defining the variables for the derivative of free energy density function
geta1 = zeros(Nx,Ny);
geta2 = zeros(Nx,Ny);
geta3 = zeros(Nx,Ny);
geta4 = zeros(Nx,Ny);
geta5 = zeros(Nx,Ny);
geta6 = zeros(Nx,Ny);
geta7 = zeros(Nx,Ny);
geta8 = zeros(Nx,Ny);
geta9 = zeros(Nx,Ny);
geta10 = zeros(Nx,Ny);
%--------------------------------------------------------------------------------------------------
% Assigning the order parameters value for ten different grain orientation  
for i = 1:Nx
  for j = 1:Ny
      if (phi(i,j) == 1)
          eta1(i,j) = 1;
      end
      if (phi(i,j) == 2)
          eta2(i,j) = 1;
      end
      if (phi(i,j) == 3)
          eta3(i,j) = 1;
      end   
      if (phi(i,j) == 4)
          eta4(i,j) = 1;
      end
      if (phi(i,j) == 5)
          eta5(i,j) = 1;
      end
      if (phi(i,j) == 6)
          eta6(i,j) = 1;
      end   
      if (phi(i,j) == 7)
          eta7(i,j) = 1;
      end
      if (phi(i,j) == 8)
          eta8(i,j) = 1;
      end
      if (phi(i,j) == 9)
          eta9(i,j) = 1;
      end   
      if (phi(i,j) == 10)
          eta10(i,j) = 1;
      end 
  end
end
b = zeros(Nx,Ny);
for i = 1:Nx
  for j = 1:Ny
      b(i,j) = eta1(i,j)*(eta2(i,j) + eta3(i,j) + eta4(i,j) + eta5(i,j) + eta6(i,j) + eta7(i,j) + eta8(i,j) + eta9(i,j) + eta10(i,j)) ...
               + eta2(i,j)*(eta3(i,j) + eta4(i,j) + eta5(i,j) + eta6(i,j) + eta7(i,j) + eta8(i,j) + eta9(i,j) + eta10(i,j)) ...
               + eta3(i,j)*(eta4(i,j) + eta5(i,j) + eta6(i,j) + eta7(i,j) + eta8(i,j) + eta9(i,j) + eta10(i,j)) ...
               + eta4(i,j)*(eta5(i,j) + eta6(i,j) + eta7(i,j) + eta8(i,j) + eta9(i,j) + eta10(i,j)) ...
               + eta5(i,j)*(eta6(i,j) + eta7(i,j) + eta8(i,j) + eta9(i,j) + eta10(i,j)) ...
               + eta6(i,j)*(eta7(i,j) + eta8(i,j) + eta9(i,j) + eta10(i,j)) ...
               + eta7(i,j)*(eta8(i,j) + eta9(i,j) + eta10(i,j)) ...
               + eta8(i,j)*(eta9(i,j) + eta10(i,j)) ...
               + eta9(i,j)*eta10(i,j);
  end
end
%--------------------------------------------------------------------------------------------------
dt = 1.0;                   %The time step
delkx = 2*pi/(Nx*dx);       %condition for fourier transform
delky = 2*pi/(Ny*dy);       %condition for fourier transform
A = 1.0;                    %setting the value of A
L = 1.0;                    %Setting the value of relaxation coefficient
kappa = 1.0;                %Setting the value of kappa (considered non-dimensional here)
%--------------------------------------------------------------------------------------------------
%                                   TEMPORAL EVOLUTION LOOP
%--------------------------------------------------------------------------------------------------
for temp = start_time : end_time
  for i = 1 : Nx
      for j = 1 : Ny
            geta1(i,j) = -eta1(i,j) + eta1(i,j)*eta1(i,j)*eta1(i,j)*eta1(i,j) ... 
                        + 2*eta1(i,j)*(eta2(i,j) + eta3(i,j) + eta4(i,j) + eta5(i,j) + eta6(i,j) + eta7(i,j) + eta8(i,j) + eta9(i,j) + eta10(i,j));
            geta2(i,j) = -eta2(i,j) + eta2(i,j)*eta2(i,j)*eta2(i,j)*eta2(i,j) ... 
                        + 2*eta2(i,j)*(eta1(i,j) + eta3(i,j) + eta4(i,j) + eta5(i,j) + eta6(i,j) + eta7(i,j) + eta8(i,j) + eta9(i,j) + eta10(i,j));
            geta3(i,j) = -eta3(i,j) + eta3(i,j)*eta3(i,j)*eta3(i,j)*eta3(i,j) ... 
                        + 2*eta3(i,j)*(eta2(i,j) + eta1(i,j) + eta4(i,j) + eta5(i,j) + eta6(i,j) + eta7(i,j) + eta8(i,j) + eta9(i,j) + eta10(i,j));
            geta4(i,j) = -eta4(i,j) + eta4(i,j)*eta4(i,j)*eta4(i,j)*eta4(i,j) ... 
                        + 2*eta4(i,j)*(eta2(i,j) + eta3(i,j) + eta1(i,j) + eta5(i,j) + eta6(i,j) + eta7(i,j) + eta8(i,j) + eta9(i,j) + eta10(i,j));
            geta5(i,j) = -eta5(i,j) + eta5(i,j)*eta5(i,j)*eta5(i,j)*eta5(i,j) ... 
                        + 2*eta5(i,j)*(eta2(i,j) + eta3(i,j) + eta4(i,j) + eta1(i,j) + eta6(i,j) + eta7(i,j) + eta8(i,j) + eta9(i,j) + eta10(i,j));
            geta6(i,j) = -eta6(i,j) + eta6(i,j)*eta6(i,j)*eta6(i,j)*eta6(i,j) ... 
                        + 2*eta6(i,j)*(eta2(i,j) + eta3(i,j) + eta4(i,j) + eta5(i,j) + eta1(i,j) + eta7(i,j) + eta8(i,j) + eta9(i,j) + eta10(i,j));
            geta7(i,j) = -eta7(i,j) + eta7(i,j)*eta7(i,j)*eta7(i,j)*eta7(i,j) ... 
                        + 2*eta7(i,j)*(eta2(i,j) + eta3(i,j) + eta4(i,j) + eta5(i,j) + eta6(i,j) + eta1(i,j) + eta8(i,j) + eta9(i,j) + eta10(i,j));
            geta8(i,j) = -eta8(i,j) + eta8(i,j)*eta8(i,j)*eta8(i,j)*eta8(i,j) ... 
                        + 2*eta8(i,j)*(eta2(i,j) + eta3(i,j) + eta4(i,j) + eta5(i,j) + eta6(i,j) + eta7(i,j) + eta1(i,j) + eta9(i,j) + eta10(i,j));
            geta9(i,j) = -eta9(i,j) + eta9(i,j)*eta9(i,j)*eta9(i,j)*eta9(i,j) ... 
                        + 2*eta9(i,j)*(eta2(i,j) + eta3(i,j) + eta4(i,j) + eta5(i,j) + eta6(i,j) + eta7(i,j) + eta8(i,j) + eta1(i,j) + eta10(i,j));
            geta10(i,j) = -eta10(i,j) + eta10(i,j)*eta10(i,j)*eta10(i,j)*eta10(i,j) ... 
                        + 2*eta10(i,j)*(eta2(i,j) + eta3(i,j) + eta4(i,j) + eta5(i,j) + eta6(i,j) + eta7(i,j) + eta8(i,j) + eta9(i,j) + eta1(i,j));        
      end
  end
  eta1hat = fft2(eta1);
  eta2hat = fft2(eta2);
  eta3hat = fft2(eta3);
  eta4hat = fft2(eta4);
  eta5hat = fft2(eta5);
  eta6hat = fft2(eta6);
  eta7hat = fft2(eta7);
  eta8hat = fft2(eta8);
  eta9hat = fft2(eta9);
  eta10hat = fft2(eta10);
  geta1hat = fft2(geta1);
  geta2hat = fft2(geta2);
  geta3hat = fft2(geta3);
  geta4hat = fft2(geta4);
  geta5hat = fft2(geta5);
  geta6hat = fft2(geta6);
  geta7hat = fft2(geta7);
  geta8hat = fft2(geta8);
  geta9hat = fft2(geta9);
  geta10hat = fft2(geta10);
  %--------------------------------------------------------------------------------------------------
  %                                     EVOLUTION EQUATION
  %--------------------------------------------------------------------------------------------------
  for i = 1:Nx
      if ((i-1) <= halfNx)  % Implementing the Periodic Boundary Condition for the x direction
      kx = (i-1)*delkx;     
      else
      kx = (i-1-Nx)*delkx;
      end
        for j = 1:Ny
            if ((j-1) <= halfNy)    % Implementing the Periodic Boundary Condition for the y direction
            ky = (j-1)*delky;       
            else
            ky = (j-1-Ny)*delky;
            end
            k2 = kx*kx + ky*ky;
           
            eta1hat(i,j) = (eta1hat(i,j) - L*dt*geta1hat(i,j))/(1 + 2*L*kappa*k2*dt);
            eta2hat(i,j) = (eta2hat(i,j) - L*dt*geta2hat(i,j))/(1 + 2*L*kappa*k2*dt);
            eta3hat(i,j) = (eta3hat(i,j) - L*dt*geta3hat(i,j))/(1 + 2*L*kappa*k2*dt);
            eta4hat(i,j) = (eta4hat(i,j) - L*dt*geta4hat(i,j))/(1 + 2*L*kappa*k2*dt);
            eta5hat(i,j) = (eta5hat(i,j) - L*dt*geta5hat(i,j))/(1 + 2*L*kappa*k2*dt);
            eta6hat(i,j) = (eta6hat(i,j) - L*dt*geta6hat(i,j))/(1 + 2*L*kappa*k2*dt);
            eta7hat(i,j) = (eta7hat(i,j) - L*dt*geta7hat(i,j))/(1 + 2*L*kappa*k2*dt);
            eta8hat(i,j) = (eta8hat(i,j) - L*dt*geta8hat(i,j))/(1 + 2*L*kappa*k2*dt);
            eta9hat(i,j) = (eta9hat(i,j) - L*dt*geta9hat(i,j))/(1 + 2*L*kappa*k2*dt);
            eta10hat(i,j) = (eta10hat(i,j) - L*dt*geta10hat(i,j))/(1 + 2*L*kappa*k2*dt);
          
        end % Ending the j loop
  end   % Ending the i loop
  %--------------------------------------------------------------------------------------------------
  eta1 = real(ifft2(eta1hat));
  eta2 = real(ifft2(eta2hat));
  eta3 = real(ifft2(eta3hat));
  eta4 = real(ifft2(eta4hat));
  eta5 = real(ifft2(eta5hat));
  eta6 = real(ifft2(eta6hat));
  eta7 = real(ifft2(eta7hat));
  eta8 = real(ifft2(eta8hat));
  eta9 = real(ifft2(eta9hat));
  eta10 = real(ifft2(eta10hat));
  %--------------------------------------------------------------------------------------------------
  if (rem(temp, time_step) == 0)
       for i = 1:Nx
          for j = 1:Ny
               b(i,j) = eta1(i,j)*(eta2(i,j) + eta3(i,j) + eta4(i,j) + eta5(i,j) + eta6(i,j) + eta7(i,j) + eta8(i,j) + eta9(i,j) + eta10(i,j)) ...
               + eta2(i,j)*(eta3(i,j) + eta4(i,j) + eta5(i,j) + eta6(i,j) + eta7(i,j) + eta8(i,j) + eta9(i,j) + eta10(i,j)) ...
               + eta3(i,j)*(eta4(i,j) + eta5(i,j) + eta6(i,j) + eta7(i,j) + eta8(i,j) + eta9(i,j) + eta10(i,j)) ...
               + eta4(i,j)*(eta5(i,j) + eta6(i,j) + eta7(i,j) + eta8(i,j) + eta9(i,j) + eta10(i,j)) ...
               + eta5(i,j)*(eta6(i,j) + eta7(i,j) + eta8(i,j) + eta9(i,j) + eta10(i,j)) ...
               + eta6(i,j)*(eta7(i,j) + eta8(i,j) + eta9(i,j) + eta10(i,j)) ...
               + eta7(i,j)*(eta8(i,j) + eta9(i,j) + eta10(i,j)) ...
               + eta8(i,j)*(eta9(i,j) + eta10(i,j)) ...
               + eta9(i,j)*eta10(i,j);
          end
       end
      mesh(b);
      colorbar;
      colormap("jet");
      title(sprintf("time %d", temp));
      xlim([1 Nx]);
      ylim([1 Ny]);
      xticks([]);
      yticks([]);
      xticklabels([]);
      yticklabels([]);   
      view(2); 
      box on;
      ax = gca;
      ax.LineWidth = 2;
      set(gcf, 'PaperUnits', 'inches', 'PaperPosition', [0 0 4 3]);
      OutputFileName = sprintf("image%d.png", temp);
      print(OutputFileName, '-dpng', '-r256');
  end
end
toc
disp('The code execution has finished');
%--------------------------------------------------------------------------------------------------
%                                       END OF CODE
%--------------------------------------------------------------------------------------------------
