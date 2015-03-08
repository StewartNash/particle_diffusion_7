%PLEASE SUBMIT MODIFICATIONS AND IMPROVEMENTS!

%File: ringfield.m (MATLAB)
%Version: 0.0
%Author: Stewart Nash
%Date: March 7, 2015
%Description: 2-Dimensional simulation of non-interacting particles originating in region with constant velocity. A reflecting boundary is present. A deflecting ring is present in the center.


%Note: Particle index means a consecutive list of integers starting with 1.
%This function is an arbitrary resonator field or coupling given by g(r,x)=g0*exp(-ar)*exp(-ikx) where x = theta and k = k0 * r0.

function output_field = ringfield(input_radius, input_angle, input_field, input_wavelength, input_wavevector, input_position, input_center)

	C_PI =  3.14159265358979323846264338327950288419716939937510582;

	alpha = 2 * C_PI / input_wavelength;
	k = input_wavevector;
	x0 = input_position(1) - input_center(1);
	y0 = input_position(2) - input_center(2);
	r = sqrt(x0^2 + y0^2);
	m0 = y0 / x0;
	
	while (input_angle > 2 * C_PI)
		input_angle = input_angle - 2 * C_PI;
	end
	while (input_angle < 0)
		input_angle = input_angle + 2 * C_PI;
	end	
	
	if (y0 >= 0 && x0 >= 0) %Quadrant I in Cartesian space
		theta = atan(m0);
	elseif (y0 >= 0 && x0 <= 0) %Quadrant II in Cartesian space
		theta = atan(m0) + C_PI;
	elseif (y0 <= 0 && x0 <= 0) %Quadrant III in Cartesian space
		theta = atan(m0) + C_PI;
	elseif (y0 <= 0 && x0 >= 0) %Quadrant IV in Cartesian space
		theta = atan(m0) + 2 * C_PI;
	else
		fprintf('ERROR (ringfield.m): Error in angle determination. Results may not be accurate.\n');
		theta = atan(m0);
	end
	
	theta = theta + input_angle;
	
	while (theta > 2 * C_PI)
		theta = theta - 2 * C_PI;
	end
	while (theta < 0)
		theta = theta + 2 * C_PI;
	end	
	
	x = input_radius * theta;		
	
	if r > input_radius
		r = r - input_radius;
	else
		r = 0;
	end
	
	output_field = input_field * exp(- alpha * r);
	output_field = output_field * exp(- 1i * k * x); 
end