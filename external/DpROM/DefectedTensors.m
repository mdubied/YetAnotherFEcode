% tensor corrections, ONLY for defect-formulations

function [Q2, Q3, Q4, Q3t, Q4t] = DefectedTensors(tensors, xi)

nd = length(xi);

% tensors computed on the nominal volume
Q2 = tensors.Q2{1};
Q3 = tensors.Q3{1};
Q4 = tensors.Q4{1};
Q3d = tensors.Q3d{1};
Q4d = tensors.Q4d{1};
Q5d = tensors.Q5d{1};
Q4dd = tensors.Q4dd{1};
Q5dd = tensors.Q5dd{1};
Q6dd = tensors.Q6dd{1};

% apply volume correction to integrate over the defected volume
if tensors.volume == 1
    for dd = 2 : nd+1
        Q2   = Q2   + tensors.Q2{dd}   * xi(dd-1);
        Q3   = Q3   + tensors.Q3{dd}   * xi(dd-1);
        Q4   = Q4   + tensors.Q4{dd}   * xi(dd-1);
        Q3d  = Q3d  + tensors.Q3d{dd}  * xi(dd-1);
        Q4d  = Q4d  + tensors.Q4d{dd}  * xi(dd-1);
        Q5d  = Q5d  + tensors.Q5d{dd}  * xi(dd-1);
        Q4dd = Q4dd + tensors.Q4dd{dd} * xi(dd-1);
        Q5dd = Q5dd + tensors.Q5dd{dd} * xi(dd-1);
        Q6dd = Q6dd + tensors.Q6dd{dd} * xi(dd-1);
    end
end

if isscalar(xi)
    % slightly different syntax if xi is a scalar
    Q3d  = Q3d*xi;
    Q4d  = ttv(Q4d,xi,3);
    Q5d  = ttv(Q5d,xi,4);
    Q4dd = Q4dd*xi^2;
    Q5dd = ttv(Q5dd*xi,xi,3);
    Q6dd = ttv(ttv(Q6dd,xi,5),xi,3);
else
    Q3d  = ttv(Q3d,xi,3);
    Q4d  = ttv(Q4d,xi,3);
    Q5d  = ttv(Q5d,xi,4);
    Q4dd = ttv(ttv(Q4dd,xi,4),xi,3);
    Q5dd = ttv(ttv(Q5dd,xi,5),xi,3);
    Q6dd = ttv(ttv(Q6dd,xi,5),xi,3);
end

Q2D = double(Q3d + Q4dd);
Q3D = double(Q4d + Q5dd);
Q4D = double(Q5d + Q6dd);

Q2 = Q2 + Q2D;
Q3 = Q3 + Q3D;
Q4 = Q4 + Q4D;

% for the tangent stiffness matrix
Q3t = Q3 + permute(Q3, [1 3 2]); 
Q4t = Q4 + permute(Q4, [1 3 2 4]) + permute(Q4, [1 4 2 3]);
