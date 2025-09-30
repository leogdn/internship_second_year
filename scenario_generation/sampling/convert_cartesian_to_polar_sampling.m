function [trajectoryPolar] = convert_cartesian_to_polar_sampling(trajectoryCartesian)

[d,T] = size(trajectoryCartesian);
trajectoryPolar = zeros(T,d);

inds_coords = [1 3];
inds_velo = [2 4];

for t=1:T
    %[trajectoryPolar(i,2),trajectoryPolar(i,1)] = cart2pol(trajectoryCartesian(i,:));
    x = trajectoryCartesian(inds_coords(1),t);
    y = trajectoryCartesian(inds_coords(2),t);
    vx = trajectoryCartesian(inds_velo(1),t);
    vy = trajectoryCartesian(inds_velo(2),t);
    [trajectoryPolar(t,3), trajectoryPolar(t,1)] = cart2pol(x,y);
    dr = (x*vx+y*vy)/trajectoryPolar(t,inds_coords(2));
    dtheta = (x*vy-y*vx)/(trajectoryPolar(t,inds_coords(2))^2);
    trajectoryPolar(t,2) = dr;
    trajectoryPolar(t,4) = dtheta;
end

end
