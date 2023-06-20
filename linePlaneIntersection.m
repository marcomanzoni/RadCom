function intersectionPoint = linePlaneIntersection(la,lb,p0,p1,p2)
%LINEPLANEINTERSECTION Calculate the intersection point of a line with a
%plane. A line is defined by the equation la + (lb-la)*t where la and lb
%are 3D vectors representing points belonging to the line. t is a real scalar. 
% The plane on the other hand, is %defined by the equation p0 + (p1-p0)*u + (p2-p0)*v 
% where p0,p1,p2 are three points on the plane and u,v are real numbers.
% Equate the two equations and solve for t,u,v.


dato = la-p0;
H = [-(lb-la), p1-p0, p2-p0];

% t = temp(1);
% u = temp(2);
% v = temp(3);
temp = H\dato; 

intersectionPoint = p0 + (p1-p0)*temp(2) + (p2-p0)*temp(3);

end

