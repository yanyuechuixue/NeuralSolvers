import sys
sys.path.append('../..')  # PINNFramework etc.

import numpy
import numpy as np
import scipy.io
from pyDOE import lhs
import torch
from torch import Tensor, ones, stack, load
from torch.autograd import grad
from torch.utils.data import Dataset
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec

sys.path.append('NeuralSolvers/')  # PINNFramework etc.
import PINNFramework as pf



class InitialConditionDataset(Dataset):

    def __init__(self, n0):
        """
        Constructor of the boundary condition dataset

        Args:
          n0 (int)
        """
        super(type(self)).__init__()
        npx,npy,npz,npvalue = np.loadtxt("./my_test_points/TwoPunctures_0.5.0.5/grid_y_axis.dat", delimiter=" ", unpack=True)

        
        Exact = npvalue

        X, T = np.meshgrid(x, t)
        xx1 = np.hstack((X[0:1, :].T, T[0:1, :].T))
        uu1 = Exact[0:1, :].T
        xx2 = np.hstack((X[:, 0:1], T[:, 0:1]))
        uu2 = Exact[:, 0:1]
        xx3 = np.hstack((X[:, -1:], T[:, -1:]))
        uu3 = Exact[:, -1:]

        X_u_train = np.vstack([xx1, xx2, xx3])
        u_train = np.vstack([uu1, uu2, uu3])

        idx = np.random.choice(X_u_train.shape[0], n0, replace=False)
        self.X_u_train = 
        self.u_train = u_train[idx, :]

    def __len__(self):
        """
        There exists no batch processing. So the size is 1
        """
        return 1

    def __getitem__(self, idx):
        x = self.X_u_train
        y = self.u_train

        return Tensor(x).float(), Tensor(y).float()
    
# Domain bounds
lb = np.array([-1, 0.0])
ub = np.array([1.0, 1.0])

nu = 0.01 / np.pi
noise = 0.0

N_u = 100
N_f = 10000

# initial condition
ic_dataset = InitialConditionDataset(n0=N_u)
initial_condition = pf.InitialCondition(ic_dataset, name='Initial condition')

#sampler
sampler = pf.LHSSampler()
#sampler = pf.RandomSampler()

# geometry
geometry = pf.NDCube(lb,ub,N_f,N_f,sampler)

# define our equations here. we can ref the wave_eq below.
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


# define underlying PDE
def burger1D(x, u):
    grads = ones(u.shape, device=u.device)  # move to the same device as prediction
    grad_u = grad(u, x, create_graph=True, grad_outputs=grads)[0]
    # calculate first order derivatives
    u_x = grad_u[:, 0]
    u_t = grad_u[:, 1]

    grads = ones(u_x.shape, device=u.device)  # move to the same device as prediction
    # calculate second order derivatives
    grad_u_x = grad(u_x, x, create_graph=True, grad_outputs=grads)[0]
    u_xx = grad_u_x[:, 0]

    # reshape for correct behavior of the optimizer
    u_x = u_x.reshape(-1, 1)
    u_t = u_t.reshape(-1, 1)
    u_xx = u_xx.reshape(-1, 1)

    f = u_t + u * u_x - (0.01 / np.pi) * u_xx
    return f

pde_loss = pf.PDELoss(geometry, burger1D, name='1D Burgers equation')
# create model
model = pf.models.MLP(input_size=2, output_size=1,
                      hidden_size=40, num_hidden=8, lb=lb, ub=ub, activation=torch.tanh)
# create PINN instance
pinn = pf.PINN(model, 2, 1, pde_loss, initial_condition, [], use_gpu=True)

logger = pf.WandbLogger("1D Burgers equation pinn", args = {})

# train pinn
pinn.fit(50000, checkpoint_path='checkpoint.pt', restart=True, logger=logger, lbfgs_finetuning=False, writing_cycle=1000)
logger.log_ptfiles()