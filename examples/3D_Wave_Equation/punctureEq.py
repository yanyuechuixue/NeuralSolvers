import torch
from torch.autograd import grad

rBD = 1 # computational domains
bhN = 2 # black hole numbers
masses = [0.128488, 0.128488] # masses of black holes
Porgs = [[0,0.2,0],[0,-0.2,0]] # original positions of black holes
Pmoms = [[-0.04,-0.001,0],[0.04,0.001,0]] # original momenta of black holes
Spins = [[0,0,0],[0,0,0]] # original spins of black holes
HLF = 0.5 
ZEO = 0.0 
ONE = 1.0 
THR = 3.0


def twopuncture_eq(x, u):

    grads = torch.ones(u.shape, device=u.device)  # move to the same device as prediction

    grad_u = grad(u, x, create_graph=True, grad_outputs=grads)[0]  # (x,y,z)

    xx = x[:, 0]
    yy = x[:, 1]
    zz = x[:, 2]

    u_x = grad_u[:, 0]
    u_y = grad_u[:, 1]
    u_z = grad_u[:, 2]

    grads = torch.ones(u_z.shape, device=u_z.device) # update for shapes

    # calculate second order derivatives
    u_zz = grad(u_z, x, create_graph=True, grad_outputs=grads)[0][:, 2]
    u_yy = grad(u_y, x, create_graph=True, grad_outputs=grads)[0][:, 1]
    u_xx = grad(u_x, x, create_graph=True, grad_outputs=grads)[0][:, 0]




    # ------------------------------------------------------------------------------
    # first calculate black hole 1
    blackholeindex = 0
    nx = xx - Porgs[blackholeindex][0]
    ny = yy - Porgs[blackholeindex][1]
    nz = zz - Porgs[blackholeindex][2]
    M = masses[blackholeindex]

    Px,Py,Pz = Pmoms[blackholeindex]
    Sx,Sy,Sz = Spins[blackholeindex]

    rr = torch.sqrt(nx**2+ny**2+nz**2) # distance to the black hole, after this line, nx,ny,nz will not be use anymore. we will use nx,ny,nz symbol to indicate the unit vector.

    nx = nx/rr # unit vector
    ny = ny/rr
    nz = nz/rr

    psi = ONE+HLF*M/rr
    tmp = Px * nx + Py * ny + Pz * nz

    Axx = ((HLF *( Px * nx + nx * Px - ( ONE - nx * nx )* tmp ) + 
        ( nx * Sy * nz - nx * Sz * ny + nx * Sy * nz - nx * Sz * ny ) / rr ) * 
        THR / ( rr * rr ))

    Ayy = ((HLF *( Py * ny + ny * Py - ( ONE - ny * ny )* tmp ) + 
        ( ny * Sz * nx - ny * Sx * nz + ny * Sz * nx - ny * Sx * nz ) / rr ) * 
        THR / ( rr * rr ))

    Azz = ((HLF *( Pz * nz + nz * Pz - ( ONE - nz * nz )* tmp ) + 
        ( nz * Sx * ny - nz * Sy * nx + nz * Sx * ny - nz * Sy * nx ) / rr ) * 
        THR / ( rr * rr ))

    Axy = ((HLF *( Px * ny + nx * Py + nx * ny * tmp ) + 
        ( nx * Sz * nx - nx * Sx * nz + ny * Sy * nz - ny * Sz * ny ) / rr ) * 
        THR / ( rr * rr ))

    Axz = ((HLF *( Px * nz + nx * Pz + nx * nz * tmp ) + 
        ( nx * Sx * ny - nx * Sy * nx + nz * Sy * nz - nz * Sz * ny ) / rr ) * 
        THR / ( rr * rr ))

    Ayz = ((HLF *( Py * nz + ny * Pz + ny * nz * tmp ) + 
        ( ny * Sx * ny - ny * Sy * nx + nz * Sz * nx - nz * Sx * nz ) / rr ) * 
    THR / ( rr * rr ))

    # blackhole 2

    blackholeindex = 1
    nx = xx - Porgs[blackholeindex][0]
    ny = yy - Porgs[blackholeindex][1]
    nz = zz - Porgs[blackholeindex][2]
    M = masses[blackholeindex]

    Px,Py,Pz = Pmoms[blackholeindex]
    Sx,Sy,Sz = Spins[blackholeindex]

    rr = torch.sqrt(nx**2+ny**2+nz**2) # distance to the black hole, after this line, nx,ny,nz will not be use anymore. we will use nx,ny,nz symbol to indicate the unit vector.

    nx = nx/rr # unit vector
    ny = ny/rr
    nz = nz/rr

    psi += HLF*M/rr
    tmp = Px * nx + Py * ny + Pz * nz

    Axx += ((HLF *( Px * nx + nx * Px - ( ONE - nx * nx )* tmp ) + 
        ( nx * Sy * nz - nx * Sz * ny + nx * Sy * nz - nx * Sz * ny ) / rr ) * 
        THR / ( rr * rr ))
    Ayy += ((HLF *( Py * ny + ny * Py - ( ONE - ny * ny )* tmp ) + 
        ( ny * Sz * nx - ny * Sx * nz + ny * Sz * nx - ny * Sx * nz ) / rr ) * 
        THR / ( rr * rr ))
    Azz += ((HLF *( Pz * nz + nz * Pz - ( ONE - nz * nz )* tmp ) + 
        ( nz * Sx * ny - nz * Sy * nx + nz * Sx * ny - nz * Sy * nx ) / rr ) * 
        THR / ( rr * rr ))
    Axy += ((HLF *( Px * ny + nx * Py + nx * ny * tmp ) + 
        ( nx * Sz * nx - nx * Sx * nz + ny * Sy * nz - ny * Sz * ny ) / rr ) * 
        THR / ( rr * rr ))
    Axz += ((HLF *( Px * nz + nx * Pz + nx * nz * tmp ) + 
        ( nx * Sx * ny - nx * Sy * nx + nz * Sy * nz - nz * Sz * ny ) / rr ) * 
        THR / ( rr * rr ))
    Ayz += ((HLF *( Py * nz + ny * Pz + ny * nz * tmp ) + 
        ( ny * Sx * ny - ny * Sy * nx + nz * Sz * nx - nz * Sx * nz ) / rr ) * 
        THR / ( rr * rr ))

    rhs = ((Axx*Axx+Ayy*Ayy+Azz*Azz+
        2*(Axy*Axy+Axz*Axz+Ayz*Ayz))*
        torch.power(psi+u, -7)/8)

    return -(u_xx + u_yy + u_zz) - rhs


if 0:
    def wave_eq(x, u):

        grads = torch.ones(u.shape, device=u.device)  # move to the same device as prediction

        grad_u = grad(u, x, create_graph=True, grad_outputs=grads)[0]  # (z, y, x, t)

        u_z = grad_u[:, 0]
        u_y = grad_u[:, 1]
        u_x = grad_u[:, 2]
        u_t = grad_u[:, 3]

        grads = torch.ones(u_z.shape, device=u_z.device) # update for shapes

        # calculate second order derivatives
        u_zz = grad(u_z, x, create_graph=True, grad_outputs=grads)[0][:, 0]  # (z, y, x, t)
        u_yy = grad(u_y, x, create_graph=True, grad_outputs=grads)[0][:, 1]
        u_xx = grad(u_x, x, create_graph=True, grad_outputs=grads)[0][:, 2]
        u_tt = grad(u_t, x, create_graph=True, grad_outputs=grads)[0][:, 3]

        f_u = u_tt - (u_zz + u_yy + u_xx)
        #relu6 = torch.nn.ReLU6()
        #propagation_error = float(1./6.) * relu6(u_y*u_t)
        return f_u
