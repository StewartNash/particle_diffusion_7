%PLEASE SUBMIT MODIFICATIONS AND IMPROVEMENTS!

%File: ringimage.m (MATLAB)
%Version: 0.0
%Author: Stewart Nash
%Date: March 7, 2015
%Description: Produce Truecolor image of domain for particles in a domain with a central ring with a field. Field is in red while particles and mask are in white.

%Note: Particle index means a consecutive list of integers starting with 1.

%>>>>Constants<<<<
%MAX_CHANNEL - Maximum channel value for Truecolor image channel
%CHANNEL_NO - Number of Truecolor image channels

%>>>>Configuration<<<<
%s_LOGSCALE - To use logarithmic field visualization, set to 1, to use linear set to 0

%>>>>Parameters<<<<
%input_values - An array of ring values as given in particle_diff_7.m

function output_image = ringimage(input_size, input_position, domain_size, input_values, input_flags)

	MAX_CHANNEL = 1.0;
	CHANNEL_NO = 3;
    
    s_LOGSCALE = 1;
	
	ring_radius = input_values(1);
	max_value = input_values(3);
    domain_center = domain_size / 2.0;
	
%	domain_image = zeros(domain_size);

	red_channel = zeros(domain_size);
	green_channel = zeros(domain_size);
	blue_channel = zeros(domain_size);

	%The following loop does not need to be run for every call of this function if the field is not time varying. Anyways, there is no variable to account for time variation in this function or the current field.
	%Please consider calculating the field only once, and saving this image channel. The present setup slows down the code considerably.	
    for j = 1 : domain_size(1)
        for k = 1 : domain_size(2)
            red_channel(j, k) = ringfield(input_values(1), input_values(2), input_values(3), input_values(4), input_values(5), [j, k], domain_center);

        end
    end
    
	%The sinusoidal component may have an imaginary component as it is given as exp(-ikx). Also, only the magnitude of the field is important, not the sign.
    red_channel = real(red_channel);
    red_channel = abs(red_channel);

	%Graph either a linear scale or a logarithmic scale.
    if s_LOGSCALE == 0
		%Scale values to the maximum allowed value.
		%This algorithm assumes that the maximum field input paramter gives the maximum value. For an alternate formulation, see the 'else' statement.
        red_channel = red_channel / max_value;
        red_channel = red_channel * MAX_CHANNEL;
    elseif s_LOGSCALE == 1
        temporary_1 = red_channel;
        red_channel = log(red_channel);
		%There is no finite logarithm of zero. Set this number to zero so we can obtain a true minimum.
        for j = 1 : domain_size(1)
            for k = 1 : domain_size(2)
                if temporary_1(j, k) == 0
                    red_channel(j, k) = 0;
                end
            end
        end
        temporary = min(red_channel);
        minimum_value = min(temporary);
		%Make the minimum value zero, eliminating negative numbers.
        red_channel = red_channel - minimum_value;
		%Set the zeros of the field to zero, the minimum value.
        for j = 1 : domain_size(1)
            for k = 1 : domain_size(2)
                if temporary_1(j, k) == 0
                    red_channel(j, k) = 0;
                end
            end
        end
		%Scale values to the maximum allowed value.
        temporary = max(red_channel);
        maximum_value = max(temporary);
        red_channel = red_channel / maximum_value;
        red_channel = red_channel * MAX_CHANNEL;
    else
        frprintf('ERROR (ringimage.m): Incorrect setting of graph scale s_LOGSCALE.');
        temporary = max(red_channel);
        maximum_value = max(temporary);
        red_channel = red_channel / maximum_value;
        red_channel = red_channel * MAX_CHANNEL;
    end
	%Clip values to maximum. However, values should not have exceeded maximum.
    for j = 1 : domain_size(1)
        for k = 1 : domain_size(2)
            if red_channel(j, k) > MAX_CHANNEL
                red_channel(j, k) = MAX_CHANNEL;
            end
        end
    end
	
	%Create mask and image of particles.
	domain_mask = updatemask(domain_size, ring_radius);	
	domain_image = updatedomainring(input_position, input_size, domain_size, input_flags);	
	domain_image = enlarge(domain_image, domain_size);
	domain_image = enlarge(domain_image, domain_size);
	domain_logical = or(domain_image, domain_mask);

	%Eliminate field visualization within the ring or disk.
	for j = 1 : domain_size(1)
		for k = 1 : domain_size(2)
			if ((j - domain_center(1))^2 + (k - domain_center(2))^2) < ring_radius^2
				red_channel(j, k) = 0;
				green_channel(j, k) = 0;
				blue_channel(j, k) = 0;
			end
		end
	end	
	%Set color to white for set pixels in mask.
	for j = 1 : domain_size(1)
		for k = 1 : domain_size(2)
			if domain_logical(j, k)
				red_channel(j, k) = MAX_CHANNEL;
				green_channel(j, k) = MAX_CHANNEL;
				blue_channel(j, k) = MAX_CHANNEL;
			end
		end
    end
    %Concatenate channels to create a Truecolor image.
	output_image = cat(CHANNEL_NO, red_channel, green_channel, blue_channel);
    
end