import numpy as np
import matplotlib.pyplot as plt

from scipy.constants import physical_constants, alpha

au_fs = physical_constants["atomic unit of time"][0] * 1.e+15
au_ev = physical_constants["atomic unit of energy"][0] / physical_constants["electron volt"][0]
au_aa = physical_constants["atomic unit of length"][0] * 1.e+10
au_va = physical_constants["atomic unit of electric field"][0] * 1.e-10
au_nm = au_aa * 1.e-1
cspeed_au = 1 / alpha

from numpy import sin, cos, exp, pi, e, sqrt

from scipy import interpolate
from scipy import signal

import sys



class DataFrame:
    
    def __init__(self, name, nx, ny, nt):
        self.name = name
        self.fh = open(name, "r")
        self.nx = nx
        self.ny = ny
        self.nz = 1
        self.nf = nt
    
    def read_field_iy(self, iy, it):
        self.fh.seek(8 * 3 * self.nx * (iy + self.ny * self.nz * it))
        temp = np.fromfile(self.fh, dtype=np.float64, count=3 * self.nx)
        return temp.reshape([self.nx, 3])
    
    def read_field(self, it):
        self.fh.seek(8 * 3 * self.nx * self.ny * it)
        temp = np.fromfile(self.fh, dtype=np.float64, count=3*self.nx*self.ny*self.nz)
        return temp.reshape([self.ny, self.nx, 3])
    




class FDTDData:
    
    
    def __init__(self, inputfile):
        data = {}
        with open(inputfile) as fh:
            for line in fh:
                temp = line.split("=")
                if len(temp) == 2:
                    data[temp[0].strip()] = temp[1].strip()
        self.sysname = eval(data["sysname"])
        self.nx1_m = int(data['nx1_m'])
        self.ny1_m = int(data['ny1_m'])
        self.nz1_m = 1#int(data['nz1_m'])
        self.nx2_m = int(data['nx2_m'])
        self.ny2_m = int(data['ny2_m'])
        self.nz2_m = 1#int(data['nz2_m'])
        self.nx_m = self.nx2_m - self.nx1_m + 1
        self.ny_m = self.ny2_m - self.ny1_m + 1
        self.nz_m = self.nz2_m - self.nz1_m + 1
        self.nt = int(data['nt'])
        self.angle = float(data['angle'].replace("d", "e"))
        self.hx_m = float(data['hx_m'].replace("d", "e"))
        self.hy_m = 250.0#float(data['hy_m'].replace("d", "e"))
        self.hz_m = 250.0#float(data['hz_m'].replace("d", "e"))
        self.dt = float(data['dt'])
        self.df = DataFrame("%s_ac.bin" % self.sysname, self.nx_m, self.ny_m, self.nt)
        self.tmax = self.dt * self.nt
        self.theta = pi * self.angle / 180.0
        
        
    
    
    def var_dump(self):
        print('# self.sysname = %s' % self.sysname)
        print('# self.nx_m = %d' % self.nx_m)
        print('# self.ny_m = %d' % self.ny_m)
        print('# self.nz_m = %d' % self.nz_m)
        print('# self.nt = %d' % self.nt)
        print('# self.angle = %e' % self.angle)
        print('# self.hx_m = %e' % self.hx_m)
        print('# self.hy_m = %e' % self.hy_m)
        print('# self.hz_m = %e' % self.hz_m)
        print('# self.dt = %e' % self.dt)
        print('# self.tmax = %e' % self.tmax)
        print('# self.theta = %e' % self.theta)
    
    
    
    def interp_field_iy(self, iy, t):
        
        def shift(fs, idelta):
            ix = np.arange(len(fs), dtype=float)
            return np.interp(ix, ix + idelta, fs)
        
        if (0 <= t < self.tmax):
            it = t / self.dt
            it1, it2 = int(it), int(it) + 1
            wt1, wt2 = it2 - it, it - it1
            ft1 = self.df.read_field_iy(iy, it1)
            ft2 = self.df.read_field_iy(iy, it2)
            return wt1 * ft1 + wt2 * ft2
            
        if (t < 0):
            f = self.df.read_field_iy(iy, 0)
            t_shift = t
        else:
            f = self.df.read_field_iy(iy, self.nt)
            t_shift = t - self.tmax

        g = np.zeros_like(f)
        ix_shift = cspeed_au * t_shift / (cos(self.theta) * self.hx_m)
        for ie in [0, 1, 2]:
            g[:, ie] = shift(f[:, ie], ix_shift)
        return g
    
    def estimate_field_iy(self, iy, t, iy_piv):
        t_t = t + self.hy_m * (iy_piv - iy) * sin(fd.theta) / cspeed_au
        return self.interp_field_iy(iy_piv, t_t)
            
    def estimate_field_iy_gauss(self, iy, t, jy1, jy2, sigma):
        f_sum = 0.0
        w_sum = 0.0
        jy_c = self.ny_m * 0.5
        for jy in range(jy1, jy2):
            w = exp(-0.5 * ((jy - jy_c) / sigma) ** 2)
            w_sum += w
            f_sum += w * self.estimate_field_iy(iy, t, jy)
        return f_sum / w_sum
    
            
            
    
    
    
def export_boundary(sysname, fd, iy1, iy2, sigma):
    
    with open("%s_bc_btm.bin" % sysname, "wb") as fh:
        for it in range(1, fd.nt + 2):
            if (it % 200 == 0):
                print(it)
            fd.estimate_field_iy_gauss(-1, it * fd.dt, iy1, iy2, sigma).tofile(fh)
    
    with open("%s_bc_top.bin" % sysname, "wb") as fh:
        for it in range(1, fd.nt + 2):
            if (it % 200 == 0):
                print(it)
            fd.estimate_field_iy_gauss(fd.ny_m, it * fd.dt, iy1, iy2, sigma).tofile(fh)





fd= FDTDData("%s.inp" % sys.argv[1])
fd.var_dump()
export_boundary(sys.argv[2], fd, 8, 24, 4)
fd = FDTDData("example2.inp")
target = range(0, 6000, 1000)
fig, axis = plt.subplots(len(target), sharex=True)
for a, it in zip(axis, target):
    a.imshow(fd.df.read_field(it)[:,:,2], aspect=1, clim=[-1,1], cmap="jet", origin="lower")
    a.plot([1000]*2, [0, 32], "--k")
plt.xlim(800,1050)
plt.xlabel("Penetration Depth")
plt.savefig(sys.argv[3])
print(sys.argv[3])






