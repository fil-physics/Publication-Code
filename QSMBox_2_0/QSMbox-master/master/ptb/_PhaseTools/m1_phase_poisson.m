function [diff_poisson,poisson_left,poisson_right] = m1_phase_poisson(phase,chi,res)

grad_phase              = dgrad(phase,res);
poisson_left            = ddiv(grad_phase,res);

grad_chi                = dgrad(chi,res);
lapl_chi                = ddiv(grad_chi,res);
gradz_chi               = grad_chi;
gradz_chi(:,:,:,1:2)    = 0;
% der2z_chi               = ddiv(gradz_chi,res);
der2z_chi               = dgrad(gradz_chi(:,:,:,3),res);
der2z_chi               = der2z_chi(:,:,:,3);
poisson_right           = ((1/3)*lapl_chi-der2z_chi);
diff_poisson = poisson_left-poisson_right;

function Gx = dgrad(chi,voxel_size)
    % Discrete Gradient Using Forward Differences with the Neuman Boundary Condition
    % Chambolle, An Algorithm for Total Variation Minimization and Applications, JMIV 2004
    % Pock et al., Global Solutions of Variational Models with Convex Regularization, SIIMS 2010

    Dx = [chi(2:end,:,:); chi(end,:,:)] - chi;
    Dy = [chi(:,2:end,:), chi(:,end,:)] - chi;
    Dz = cat(3, chi(:,:,2:end), chi(:,:,end)) - chi;

    Dx = Dx/voxel_size(1);
    Dy = Dy/voxel_size(2);
    Dz = Dz/voxel_size(3);

    Gx = cat(4, Dx, Dy, Dz);
end


function div = ddiv(Gx,voxel_size)
    % Discrete Divergence Using Backward Difference with the Dirichlet Boundary Condition
    % Chambolle, An Algorithm for Total Variation Minimization and Applications, JMIV 2004

    Gx_x = Gx(:,:,:,1);
    Gx_y = Gx(:,:,:,2);
    Gx_z = Gx(:,:,:,3);

    [Mx, My, Mz] = size(Gx_x);

    Dx = [Gx_x(1:end-1,:,:); zeros(1,My,Mz)]...
        - [zeros(1,My,Mz); Gx_x(1:end-1,:,:)];

    Dy = [Gx_y(:,1:end-1,:), zeros(Mx,1,Mz)]...
        - [zeros(Mx,1,Mz), Gx_y(:,1:end-1,:)];

    Dz = cat(3, Gx_z(:,:,1:end-1), zeros(Mx,My,1))...
        - cat(3, zeros(Mx,My,1), Gx_z(:,:,1:end-1));

    Dx = Dx/voxel_size(1);
    Dy = Dy/voxel_size(2);
    Dz = Dz/voxel_size(3);

    div = -( Dx + Dy + Dz );
end

end