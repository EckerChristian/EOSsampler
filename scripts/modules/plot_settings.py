import matplotlib.pyplot as plt
import os
from modules import utils as u
from matplotlib import cm
import matplotlib.colors as mcolors
import numpy as np
from mpl_toolkits.axes_grid1 import make_axes_locatable



# Settings for different figure types
def dn_npt_settings(ax):
    ax.set_xscale('linear')
    ax.set_yscale('linear') 
    ax.minorticks_on()
    ax.set_ylim(u.DNMIN,u.DNMAX)
    ax.set_xlim(u.NPTMIN,u.NPTMAX_TWINS)
    #ax.set_xticks([1,2,3,4,5])
    #ax.set_yticks([1,2,3,4,5])
    ax.set_xticks([2,4,6,8,10])
    
    ax.set_xlabel(r'$n_{\mathrm{PT}}$ $[n_\mathrm{s}]$ ')
    #ax.set_xlabel(r'$n_{\mathrm{destab}}$ $[n_\mathrm{s}]$ ')
    ax.set_ylabel(r'$\Delta n$ $[n_\mathrm{s}]$ ')
    #ax.set_ylabel(r'$\Delta n_{\mathrm{destab}}$ $[n_\mathrm{s}]$ ')


def dngap_n1_settings(ax):
    ax.set_xscale('linear')
    ax.set_yscale('linear') 
    ax.minorticks_on()
    ax.set_ylim(u.DNGAP_MIN,u.DNGAP_MAX)
    ax.set_xlim(u.NPTMIN,u.NPTMAX_TWINS)
    #ax.set_xticks([1,2,3,4,5])
    ax.set_yticks([2,4,6,8,10])
    ax.set_xticks([2,4,6,8,10])
    ax.set_xlabel(r'$n_1$ $[n_\mathrm{s}]$ ')
    ax.set_ylabel(r'$\delta n$ $[n_\mathrm{s}]$ ')


def EOS_settings(ax):
    ax.set_xscale('log')
    ax.set_yscale('log') 
    ax.minorticks_on() 
    ax.set_xlim(u.EMIN,u.EMAX)
    ax.set_ylim(u.PMIN,u.PMAX)
    ax.set_xlabel(r'$e$ $[\mathrm{MeV/fm}^{3}]$')
    ax.set_ylabel(r'$p$ $[\mathrm{MeV/fm}^{3}]$')

def EOS_settings2(ax):
    ax.set_xscale('log')
    ax.set_yscale('log') 
    ax.minorticks_on() 
    ax.set_xlim(u.NMIN,u.NMAX)
    ax.set_ylim(u.PMIN,u.PMAX)
    ax.set_xlabel(r'$n$ $[n_\mathrm{s}]$')
    ax.set_ylabel(r'$p$ $[\mathrm{MeV/fm}^{3}]$')

def MR_settings(ax):
    ax.set_xscale('linear')
    ax.set_yscale('linear')
    ax.minorticks_on() 
    ax.set_xlim(u.RMIN, u.RMAX)
    ax.set_ylim(u.MMIN, u.MMAX)
    ax.set_xlabel(r'$R$ $[\mathrm{km}]$') 
    ax.set_ylabel(r'$M$ $[M_\odot]$')
    ax.set_xticks([8,10,12,14,16])
    #ax.set_yticks([1,2,3,4,5])

def CS2_settings(ax):
    ax.set_xscale('log')
    ax.set_yscale('linear')
    ax.minorticks_on()
    #ax.set_xlim(u.EMIN, u.EMAX)
    ax.set_xlim(u.NMIN, u.NMAX)
    ax.set_ylim(u.CS2MIN, u.CS2MAX) 
    #ax.set_xlabel(r'Energy density [MeV/fm$^{3}$]')
    ax.set_xlabel(r'$n$ $[n_\mathrm{s}]$')
    ax.set_ylabel(r'$c_\mathrm{s}^2$')
    ax.set_yticks([0, 0.2, 0.4, 0.6, 0.8])


def delta_settings(ax):
    ax.set_xscale('log')
    ax.set_yscale('linear')
    ax.minorticks_on()
    #ax.set_xlim(u.EMIN, u.EMAX)
    ax.set_xlim(u.NMIN, u.NMAX)
    ax.set_ylim(u.DMIN, u.DMAX)
    #ax.set_xlabel(r'Energy density [MeV/fm$^{3}$]')
    ax.set_xlabel(r'$n$ $[n_\mathrm{s}]$')
    ax.set_ylabel(r'$\Delta$', labelpad = -5)
    ax.set_yticks([-0.2, 0, 0.2])

def dc_settings(ax):
    ax.set_xscale('log')
    ax.set_yscale('linear')
    ax.minorticks_on()
    #ax.set_xlim(u.EMIN, u.EMAX)
    ax.set_xlim(u.NMIN, u.NMAX)
    ax.set_ylim(u.DCMIN, u.DCMAX)
    #ax.set_xlabel(r'Energy density [MeV/fm$^{3}$]')
    ax.set_xlabel(r'$n$ $[n_\mathrm{s}]$')
    ax.set_ylabel(r'$d_\mathrm{c}$')
    ax.set_yticks([0, 0.2, 0.4, 0.6, 0.8])

def dn_dngap_settings(ax): ######################################
    ax.set_xscale('linear')
    ax.set_yscale('linear')
    ax.minorticks_on() 
    ax.set_xlim(u.DNGAP_MIN, u.DNGAP_MAX)
    ax.set_ylim(u.DNMIN, u.DNMAX)
    ax.set_xlabel(r'$\delta n_\mathrm{gap}$ $[n_s]$') 
    ax.set_ylabel(r'$\Delta n$ $[n_s]$')
    ax.set_xticks([0,2,4,6,8,10])
    ax.set_yticks([2,4,6,8,10])

def nPT_n1_settings(ax):
    ax.set_xscale('linear')
    ax.set_yscale('linear')
    ax.minorticks_on() 
    ax.set_xlim(u.N1MIN, u.N1MAX)
    ax.set_ylim(u.N1MIN, u.NPTMAX)
    ax.set_xlabel(r'$n_1$ $[n_s]$') 
    ax.set_ylabel(r'$n_\mathrm{PT}$ $[n_s]$')
    ax.set_xticks([0,2,4,6,8,10])
    ax.set_yticks([5,10,15,20])

# Other plot settings
def configure_plot_style():
    """
    Make all plots look the same
    """
    plt.rcParams.update({
        "text.usetex": False,
        "font.family": "serif",
        "mathtext.fontset": "cm",
        "xtick.labelsize": 11,
        "ytick.labelsize": 11,
        "axes.labelsize": 12,
        "axes.titlesize": 12,
        "legend.fontsize": 8,
        "legend.title_fontsize": 8,
        "figure.titlesize": 16,
    })



def colormap(fig, axs, fraction=0.02, pad = 0.08):
    """
    Colormap settings for the plots.
    """
    #cmap = cm.get_cmap('plasma_r')
    #cmap = (mcolors.ListedColormap(plt.cm.Blues(np.linspace(0.3,0.95,10))))#.with_extremes(under='lightgrey'))
    #cmap = (mcolors.ListedColormap(plt.cm.plasma_r(np.linspace(0,0.85,10))).with_extremes(under='lightgrey'))
    
    colors = ['lightgrey'] + list(plt.cm.Blues(np.linspace(0.3,0.95,11)))
    cmap = mcolors.ListedColormap(colors)

    ##bounds = [0.01, 0.12, 0.19, 0.28, 0.41, 0.55, 0.69, 0.81, 0.91, 0.98, 1]
    ##bounds = [1e-4, 1e-3, 1e-2, 0.1, 0.2, 0.3, 0.4, 0.5, 1]
    #bounds = [0, 1e-4, 1e-3, 1e-2, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 1]
    bounds = [0, 1e-3, 1e-2, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.8, 1]
    norm = mcolors.BoundaryNorm(bounds, cmap.N)#, extend='min')
    cbar = fig.colorbar(plt.cm.ScalarMappable(norm=norm, cmap=cmap), orientation='vertical',  ax=axs, fraction=fraction, pad=pad, label='Normalized likelihood')
    cbar.set_ticks(bounds)
    #cbar.set_ticklabels(["0", r"$10^{-4}$", r"$10^{-3}$", r"$10^{-2}$", "0.1", "0.2", "0.3", "0.4", "0.5", "0.6", "1"])
    cbar.set_ticklabels(["0", r"$10^{-3}$",r"$10^{-2}$", "0.1", "0.2", "0.3", "0.4", "0.5", "0.6", "0.8", "1"])

    cbar.ax.tick_params(labelsize=11)

    return cmap, bounds, norm, cbar



def colormap3(fig, axs, fraction=0.02, pad = 0.08):
    """
    Colormap settings without the grey from 0 to lower bound.
    """
    cmap = (mcolors.ListedColormap(plt.cm.Blues(np.linspace(0.3,0.95,10))))#.with_extremes(under='lightgrey'))
    #bounds = [1e-4, 1e-3, 1e-2, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 1]
    bounds = [1e-2, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.8, 1]
    norm = mcolors.BoundaryNorm(bounds, cmap.N)#, extend='min')

    cbar = fig.colorbar(plt.cm.ScalarMappable(norm=norm, cmap=cmap), orientation='horizontal', ax=axs, fraction=fraction, pad=pad, label='Normalized likelihood', location='top')
    cbar.set_ticks(bounds)
    #cbar.set_ticklabels([r"$10^{-4}$", r"$10^{-3}$", r"$10^{-2}$", "0.1", "0.2", "0.3", "0.4", "0.5", "0.6", "1"])
    cbar.set_ticklabels([r"$10^{-2}$", "0.1", "0.2", "0.3", "0.4", "0.5", "0.6", "0.8", "1"])
    
    
    cbar.ax.tick_params(labelsize=11)


    return cmap, bounds, norm, cbar





def colormap_dc():
    #cmap = plt.get_cmap('RdYlBu')
    #colors = cmap([0.2, 0.75, 0.85, 1])
    cmap = (mcolors.ListedColormap(plt.cm.Blues(np.linspace(0.3,0.95,10))))
    bounds = [1e-2, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.8, 1]
    norm = mcolors.BoundaryNorm(bounds, cmap.N)#, extend='min')


    return cmap, bounds, norm




def subplot_settings(nrow_fig, ncol_fig, figsize=(5,4), wspace=0.1, hspace=0.1, sharex=True, sharey=True, squeeze=False):
    """
    Settings for subplots in a figure. Intended for figures showing the astro constraint progression.
    """

    fig, axs = plt.subplots(
        nrows=nrow_fig, ncols=ncol_fig, figsize=figsize,
        sharex=sharex, sharey=sharey, squeeze=squeeze,
    )

    plt.subplots_adjust(wspace=wspace, hspace=hspace) 


    for j in range(ncol_fig):
        for i in range(nrow_fig):
            axs[i, j].tick_params(which='major', direction='out',top="true",right="true")
            axs[i, j].tick_params(which='major', direction='in', left="true")
            axs[i, j].tick_params(which='minor', direction='in',top="true",right="true", left="true", bottom="true")

            if j > 0:   # Hide ticks on non-left columns
                axs[i, j].tick_params(axis='y', which='both', left=True, right=True, labelleft=False)

    return fig, axs

    


# Save plots
def save_plots(save=False, show=True, save_dir=".", filename="plot.png"):
    """
    Save the plot to a file.
    """
    #plt.tight_layout()

    if save:
        os.makedirs(save_dir, exist_ok=True)
        plt.savefig(os.path.join(save_dir, filename), dpi = 300, bbox_inches='tight')
        print(f"Plot saved to {os.path.join(save_dir, filename)}")
    if show:
        plt.show()
    else:
        plt.close()
