function volume = tetrahedron_volume(sep12, sep13, sep14)
%TETRAHEDRON_VOLUME Calculates the volume of a tetrahedron
%   Inputs are 3D vectors defining the radius between point 1 and the other
%   points on the tetrahedron
volume = (1/6) * abs(dot(cross(sep12, sep13),sep14));

end