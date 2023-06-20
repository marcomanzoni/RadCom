function l = pathLengthUnderSnow(Sx, Sy, Sz, px, py, pz, installation_height)
%PATHLENGTHUNDERSNOW Calculate the approximate path length under the snow.
%Calculate the intersection point of a line with a
%plane. A line is defined by the equation la + (lb-la)*t where la and lb
%are 3D vectors representing points belonging to the line. t is a real scalar. 
% The plane on the other hand, is %defined by the equation p0 + (p1-p0)*u + (p2-p0)*v 
% where p0,p1,p2 are three points on the plane and u,v are real numbers.
% Equate the two equations and solve for t,u,v.


positionSensor = [Sx; Sy; Sz];
positionTarget = [px; py; pz];
p0 = [0;0;-installation_height];
p1 = [1;0;-installation_height];
p2 = [0;1;-installation_height];

d = [positionTarget-p0];
H = [-(positionSensor-positionTarget), p1-p0, p2-p0];

temp = H\d;

% Point in space where the EM ray hits the ground. 
intersectionPoint = p0 + (p1-p0)*temp(2) + (p2-p0)*temp(3);

% Length of the path under the snow;
l = sqrt((intersectionPoint(1)-px)^2+(intersectionPoint(2)-py)^2+(intersectionPoint(3)-pz)^2);


end

