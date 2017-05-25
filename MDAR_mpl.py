from MDAR import mdar
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.widgets import Slider, Button, RadioButtons
from scipy.interpolate import interp1d

hlittle = 0.7
rhoc    = 2.775e2*hlittle**2
G       = 4.299e-6
Delta   = 200.0

#fig, ax = plt.subplots()
fig = plt.figure(1, figsize=(7,7))
plt.subplots_adjust(left=0.25, bottom=0.25)
t = np.arange(0.0, 1.0, 0.001)
a0 = 5
f0 = 3
s = a0*np.sin(2*np.pi*f0*t)
#l, = plt.plot(t, s, lw=2, color='red')
#plt.axis([0, 1, -10, 10])

ax1 = fig.add_axes([0.1,0.67,0.3,0.3])#
ax2 = fig.add_axes([0.6,0.67,0.3,0.3])#
ax3 = fig.add_axes([0.1,0.3,0.3,0.3])#
ax4 = fig.add_axes([0.6,0.3,0.3,0.3])#
#ax2 = fig.add_axes([0.4,0.35,0.65,0.2])
#ax3 = fig.add_axes([0.45,0.30,0.65,0.2])
#ax4 = fig.add_axes([0.25,0.30,0.65,0.2])

axcolor = 'lightgoldenrodyellow'
ax_Mgal      = plt.axes([0.1, 0.02, 0.8, 0.03], facecolor=axcolor)
ax_Rh        = plt.axes([0.1, 0.06, 0.8, 0.03], facecolor=axcolor)
ax_m200      = plt.axes([0.1, 0.10, 0.8, 0.03], facecolor=axcolor)
ax_fmbulge   = plt.axes([0.1, 0.14, 0.8, 0.03], facecolor=axcolor)
ax_fabulge   = plt.axes([0.1, 0.18, 0.8, 0.03], facecolor=axcolor)

ax1.set_xlabel(r'$\rm log_{10} \left ( M_{200} / M_{\odot} \right )$')
ax1.set_ylabel(r'$\rm log_{10} \left ( M_{str} / M_{\odot} \right )$')

ax2.set_xlabel(r'$\rm log_{10} \left ( r_{50} / kpc \right )$')
ax2.set_ylabel(r'$\rm log_{10} \left ( M_{str} / M_{\odot} \right )$')

ax3.set_xlabel(r'$\rm log_{10} \left ( g_{tot} / m \ s^{-2} \right )$')
ax3.set_ylabel(r'$\rm log_{10} \left ( g_{bar} / m \ s^{-2} \right )$')

ax4.set_xlabel(r'$\rm log_{10} \left ( R / \ kpc \right )$')
ax4.set_ylabel(r'$\rm log_{10} \left ( V_c / km \ s^{-1} \right )$')

mymdar = mdar()

M200 = mymdar.M200
Mgal = mymdar.Mgal
Rh   = mymdar.Rh
f_mbulge = mymdar.f_mbulge
f_abulge = mymdar.f_abulge

mymdar.plot_mdar(ax=ax3)

gdm    = mymdar.get_halo_acc(mymdar.r)
gdisk  = mymdar.get_disk_acc(mymdar.r)

plot3, = ax3.plot(np.log10(gdisk), np.log10(gdm+gdisk))

ax3.set_xlim(-13,-7)
ax3.set_ylim(-13,-7)

vcirc_bulge = mymdar.get_vcirc_bulge(mymdar.r)
vcirc_disk  = mymdar.get_vcirc_disk(mymdar.r)
vcirc_halo  = mymdar.get_vcirc_halo(mymdar.r)
vcirc_tot   = np.sqrt(vcirc_bulge**2+vcirc_disk**2+vcirc_halo**2)

plot4_disk,  = ax4.plot(mymdar.r, vcirc_disk,  label='disk')
plot4_halo,  = ax4.plot(mymdar.r, vcirc_halo,  label='halo')
plot4_bulge, = ax4.plot(mymdar.r, vcirc_bulge, label='bulge')
plot4_tot,   = ax4.plot(mymdar.r, vcirc_tot,   label='total')

ax4.set_xlim(0,100)
ax4.set_ylim(0,400)

#a.plot_galaxy(ax=ax3)


#mymdar.plot_vcirc(ax=ax4)
#ax4.legend(loc='upper left', ncol=2, frameon=False)

s_Mgal      = Slider(ax_Mgal, 'Mbar', 3, 15.0, valinit=np.log10(Mgal))
s_Rh        = Slider(ax_Rh, 'R50', 0.01, 15.0, valinit=Rh)
s_m200      = Slider(ax_m200, 'M200', 10, 14, valinit=np.log10(M200))
s_fmbulge   = Slider(ax_fmbulge, 'Bul_frac', 0.01, 1.0, valinit=f_mbulge)
s_fabulge   = Slider(ax_fabulge, 'Bul_size', 0.01, 3.0, valinit=f_abulge)


data_behroozi = np.loadtxt('m200-mstr.txt')
ax1.plot(np.log10(data_behroozi[:,0]),
         np.log10(data_behroozi[:,1]*data_behroozi[:,0]), ls='--')

galaxy_mstr_mhalo, = ax1.plot(np.log10(mymdar.M200),
                              np.log10(mymdar.Mgal), 'o')

mass_concentration = np.loadtxt('mass_concentration.txt')
get_concentration = interp1d(10**(mass_concentration[:,0]+10-np.log10(hlittle)),
                             mass_concentration[:,1])
get_galaxy_mass = interp1d(data_behroozi[:,0],
                           data_behroozi[:,1]*data_behroozi[:,0])

m200_fid = 10**np.linspace(10.2,14.0)
mgal_fid = get_galaxy_mass(m200_fid)
r200_fid = (m200_fid/(Delta*rhoc*4./3.*np.pi))**(1./3.)
c        = get_concentration(m200_fid)
rs       = r200_fid/c
rh_fid   = 0.2*rs

ax2.plot(np.log10(rh_fid), np.log10(mgal_fid), '--')
galaxy_sizes, = ax2.plot(np.log10(mymdar.Rh), np.log10(mymdar.Mgal), 'o')
        
def update(val):
    Mgal      = s_Mgal.val
    M200      = s_m200.val
    f_mbulge   = s_fmbulge.val
    f_abulge   = s_fabulge.val
    Rh        = s_Rh.val
    
    if(Mgal == 0):
        Mgal == 0
    else:
        Mgal = 10**Mgal

    if(M200 == 0):
        M200 == 0
    else:
        M200 = 10**M200
        
    mymdar.set_new_paramters(M200=M200, f_mbulge=f_mbulge, f_abulge=f_abulge, Mgal=Mgal, Rh=Rh)

    gdm    = mymdar.get_halo_acc(mymdar.r)
    gdisk  = mymdar.get_disk_acc(mymdar.r)
    plot3.set_data(np.log10(gdisk),np.log10(gdm+gdisk))

    vcirc_bulge = mymdar.get_vcirc_bulge(mymdar.r)
    vcirc_disk  = mymdar.get_vcirc_disk(mymdar.r)
    vcirc_halo  = mymdar.get_vcirc_halo(mymdar.r)
    vcirc_tot   = np.sqrt(vcirc_bulge**2+vcirc_disk**2+vcirc_halo**2)
    
    plot4_disk.set_data(mymdar.r, vcirc_disk)
    plot4_halo.set_data(mymdar.r, vcirc_halo)
    plot4_bulge.set_data(mymdar.r, vcirc_bulge)
    plot4_tot.set_data(mymdar.r, vcirc_tot)

    
    galaxy_mstr_mhalo.set_data(np.log10(M200), np.log10(Mgal))

    galaxy_sizes.set_data(np.log10(mymdar.Rh), np.log10(mymdar.Mgal))
    
    fig.canvas.draw_idle()
#sfreq.on_changed(update)
s_Mgal.on_changed(update)
s_m200.on_changed(update)
s_Rh.on_changed(update)
s_fmbulge.on_changed(update)
s_fabulge.on_changed(update)


plt.show()


