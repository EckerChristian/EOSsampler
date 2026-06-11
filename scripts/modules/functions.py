from modules import plot_settings as ps
import numpy as np
from modules import utils as u
from matplotlib.colors import LogNorm
import matplotlib as plt
import seaborn as sns
import pandas as pd

# Get plot settings 
ps.configure_plot_style()


# Auxiliary functions



def credible_intervals2(PT_data, PT_constraint, nonPT_data, nonPT_constraint, bins, interval,text):

    q = ((100 - interval) / 2) / 100
    #print(q)


    # PT data
    if PT_data is not None:
        mask_PT = PT_constraint > 0
        PT_data = PT_data[mask_PT]
        PT_constraint = PT_constraint[mask_PT]

    # nonPT data
    if nonPT_data is not None:
        mask_nonPT = nonPT_constraint > 0
        nonPT_data = nonPT_data[mask_nonPT]
        nonPT_constraint = nonPT_constraint[mask_nonPT]

    if PT_data is None:
        data = nonPT_data
        weights = nonPT_constraint
    elif nonPT_data is None:
        data = PT_data
        weights = PT_constraint
    else:
        data = np.concatenate((PT_data, nonPT_data))
        weights = np.concatenate((PT_constraint, nonPT_constraint))

    counts, bin_edges = np.histogram(data, weights=weights, bins= bins)

    lower = interp_quantile(counts, bin_edges, q)
    median = interp_quantile(counts, bin_edges, 0.5)
    upper = interp_quantile(counts, bin_edges, 1-q)

    #print(text, lower, median, upper)
    print(f"{text}_interval = [{lower}, {median}, {upper}]")


# n1, delta n, ptot, n1, delta n, ptot, family, ...  
def credible_intervals_params(n1_dataPT, dn_dataPT, constraint_dataPT, n1_dataNonPT, dn_dataNonPT, constraint_dataNonPT, bins, interval, text, family=u.EoSFamily.FAM_1, param=u.Param.n1):

    q = ((100 - interval) / 2) / 100
    #print(q)

    match family:
        case u.EoSFamily.FAM_1: # in ns
            n1_lower = 0
            n1_upper = 3
            dn_lower = 0
            dn_upper = 3
        case u.EoSFamily.FAM_2: # in ns
            n1_lower = 3
            n1_upper = 15
            dn_lower = 0
            dn_upper = 3

    match param:
        case u.Param.n1:
            PT_data = n1_dataPT
            nonPT_data = n1_dataNonPT
        case u.Param.dn:
            PT_data = dn_dataPT
            nonPT_data = dn_dataNonPT


    # PT data
    if n1_dataPT is not None and dn_dataPT is not None: # ptot > 0, n1 and delta n range
        mask_PT = (constraint_dataPT > 0) & (n1_dataPT >= n1_lower) & (n1_dataPT < n1_upper) & (dn_dataPT >= dn_lower) & (dn_dataPT < dn_upper)
        PT_data = PT_data[mask_PT]
        constraint_dataPT = constraint_dataPT[mask_PT]

    # nonPT data
    if n1_dataNonPT is not None and dn_dataNonPT is not None:
        mask_nonPT = (constraint_dataNonPT > 0) & (n1_dataNonPT >= n1_lower) & (n1_dataNonPT < n1_upper) & (dn_dataNonPT >= dn_lower) & (dn_dataNonPT < dn_upper)
        nonPT_data = nonPT_data[mask_nonPT]
        constraint_dataNonPT = constraint_dataNonPT[mask_nonPT]

    if PT_data is None:
        data = nonPT_data
        weights = constraint_dataNonPT
    elif nonPT_data is None:
        data = PT_data
        weights = constraint_dataPT
    else:
        data = np.concatenate((PT_data, nonPT_data))
        weights = np.concatenate((constraint_dataPT, constraint_dataNonPT))

    counts, bin_edges = np.histogram(data, weights=weights, bins= bins)

    lower = interp_quantile(counts, bin_edges, q)
    median = interp_quantile(counts, bin_edges, 0.5)
    upper = interp_quantile(counts, bin_edges, 1-q)

    #print(text, lower, median, upper)
    print(f"{text}_interval = [{lower}, {median}, {upper}]")



def interp_quantile(counts, bin_edges, q):
    cum = np.cumsum(counts)
    cum = cum / cum[-1] # Normalize to 1
    i = np.searchsorted(cum, q)
    
    if i == 0:
        return 0.5 * (bin_edges[0] + bin_edges[1])
    else:
        # fraction of the way through bin i
        cdf_prev = cum[i-1]
        cdf_bin = cum[i] - cdf_prev
        frac = (q - cdf_prev) / cdf_bin

        # bin edges
        left = bin_edges[i]
        right = bin_edges[i+1]

        return left + frac * (right - left)
    

def bayes_factors(likelihoods_ns, Nns, likelihoods_twins, Ntwins):

    average_likelihood_ns = likelihoods_ns / Nns
    average_likelihood_twins = likelihoods_twins / Ntwins

    bayes_factors = average_likelihood_twins / average_likelihood_ns

    print(f"bayes_factors = [{bayes_factors[0]}, {bayes_factors[1]}, {bayes_factors[2]}, {bayes_factors[3]}]")


def cuts(data, max_neos, reverse = True):
    """
    Cuts for ptot
    """
    data = sorted(data, reverse=reverse)
    data = data[:max_neos]
    
    return data, data[-1]



    

def PT_parameters(first_branch, second_branch, n_cent):
    """
    PT parameters for non-FOPT induced twins
    """

    # End of the first branch [ns]
    n1 = n_cent[first_branch][-1]
    
    # Start of the second branch [ns]
    n2 = n_cent[second_branch][0]

    # End of the second branch
    n3 = n_cent[second_branch][-1]

    # Size of the unstable part [ns]
    dn = n2 - n1

    return dn, n1, n2, n3





def find_max(M_data, stab):
    """
    Find the maximum mass of NS and Twins and stable branches. Return None for M2 and second branch for NS
    """
    first_branch = np.where(stab == 1)[0]

    if len(first_branch) < 2 or first_branch is None:
        return None, None, None, None

    first_max = first_branch[-1]
    M_1 = M_data[first_max]

    if np.max(stab) == 2:

        #first and second branch
        second_branch = np.where(stab == 2)[0]

        # If the second branch has less than 2 points, ignore it
        if len(second_branch) < 2:
            M_2 = None
            second_branch = None
        else:
            #location of the first and second max
            second_max = second_branch[-1]

            #find the first and second max
            M_2 = M_data[second_max]

        return M_1, M_2, first_branch, second_branch

    else:
        return M_1, None, first_branch, None




def find_category(stability, n1, n2, n_cent, npt):
    """
    Find the category of PT or nonPT induced twins
    """

    # There is no PT
    if npt == 0:
        return 4

    # Find the location of the last stable point
    loc_end_e = np.searchsorted(n_cent, npt) - 1 

    # nonPT twins, both branches start before PT, 2nd branch can extend beyond PT
    if n1 < npt and n2 < npt:
        return 4
    
    # PT twins, one stable point
    elif stability[loc_end_e] == 1.0 and stability[loc_end_e+1] < 0: 
        return 1
    
    # PT twins, one unstable point
    if stability[loc_end_e-1] == 1.0 and stability[loc_end_e] == -1.0 and stability[loc_end_e+1] < 0: # PT twins
        return 1
    
    # nonPT twins, destablization after PT
    elif stability[loc_end_e+1] == 1.0 and stability[loc_end_e+2] == 1.0: 
        return 2

    # nonPT twins, at least two unstable points
    elif stability[loc_end_e-1] == -1: 
        return 3

    # NonPT twins, one stable point, resemble crossover twins
    elif stability[loc_end_e] == 1 and stability[loc_end_e+1] == 1 and stability[loc_end_e+2] == -1:
        return 6
    
    # If nothing matches
    else:
        print("??")
        return 7
        


def prepare_branches(x_data, y_data, xp, fp, x_start=None, y_start=None): ##!!!!!!!!!!!!!
    """
    Prepare branches for stability plotting.
    """

    x_end, y_end = x_data[-1], y_data[-1] 

    # first branch
    if x_start is None and y_start is None: 
        # mask the range up to x_end
        mask = xp <= x_end 

    # second branch
    else: 
        # mask the range from x_start to x_end
        mask = (xp >= x_start) & (xp <= x_end) 
            
    # Logical mask for range
    xp_cut = xp[mask]
    fp_cut = fp[mask]

    # Ensure endpoints are included
    if len(xp_cut) == 0 or xp_cut[-1] < x_end:
        xp_cut = np.append(xp_cut, x_end)
        fp_cut = np.append(fp_cut, y_end)

    # Ensure starting points is included
    if x_start is not None and xp_cut[0] > x_start:
        xp_cut = np.insert(xp_cut, 0, x_start)
        fp_cut = np.insert(fp_cut, 0, y_start)


    return xp_cut, fp_cut


def astro_constraints(file, i):
    """
    Return astro constraints
    """
    params_path = f"/{i}/params"

    # CET
    pCET = file[f"{params_path}/pCET"][:][0]

    # Radio Pulsar Mass measurements
    RPM = file[f"{params_path}/pM"][:][0]  

    # GW170817
    pGW = file[f"{params_path}/pGW"][:][0]

    # Xray measurements
    Xrays = file[f"{params_path}/pXray"][:]
    #pXray = np.prod(Xrays)

    # NICER measurements from Xrays
    pNICER = Xrays[11] * Xrays[12] * Xrays[13] * Xrays[14]

    # return r, r_GW. r_GW_NICER
    return pCET, pCET*RPM, pCET*RPM*pGW, pCET*RPM*pGW*pNICER    




def plot_twins_stability(data, xp, fp, M_1, M_2, first_branch, second_branch, plot, likelihood, ax):
    """
    Plot the stability of twins
    """

    likelihood = 1  # debugging

    # find the category
    category = find_category(M_1, M_2)

    # Plot if the category is right and likeligood is nonzero
    if category == plot and likelihood > 0:

        # loop over the branches
        for where_stab in [first_branch, second_branch]:
            x_data = data[where_stab] 
            y_data = np.interp(x_data, xp, fp)

            # Second branch
            if where_stab is second_branch:
                x_start, y_start = x_data[0], y_data[0]
                xp_cut, fp_cut = prepare_branches(x_data, y_data, xp, fp, x_start=x_start, y_start=y_start)

            # First branch
            else:
                xp_cut, fp_cut = prepare_branches(x_data, y_data, xp, fp)

            # Plot
            #print("cut", xp_cut, fp_cut)
            #axs[0,0].plot(x_end, y_end, 'o', markeredgecolor='red', c='red', zorder = 2)
            #axs[0,0].plot(x_data, y_data, 'o', markeredgecolor='black', c='C'+str(i), zorder = 1)
            ax.plot(xp_cut, fp_cut, c='C'+str(i), linewidth=1,zorder = 1)#, 'o', markeredgecolor='black')
            ax.plot(xp, fp, c='grey', linewidth= 0.5, zorder = 0)



def make_grid(nx, ny, depth=None, grid=u.GridType.LINEAR, xmin=0.0, xmax=1.0, ymin=0.0, ymax=1.0):
    """
    Make linear or logarithmic grid.
    """

    match grid:
        case u.GridType.LINEAR:
            xgrid = np.linspace(xmin, xmax, nx)
            ygrid = np.linspace(ymin, ymax, ny)

            if depth is None:
                bins = np.zeros((ny-1,nx-1))
            else:
                bins = np.zeros((depth,ny-1,nx-1))

            return bins, xgrid, ygrid

        case u.GridType.LOGARITHMIC:
            if xmin <= 0 or xmax <= 0 or ymin <= 0 or ymax <= 0:
                raise ValueError("Logarithmic grid requires positive bounds")
            else:
                xgrid = [np.exp(np.log(xmin)*(1.0 - cnt/nx) + np.log(xmax)*(cnt/nx)) for cnt in range(nx)]
                ygrid = [np.exp(np.log(ymin)*(1.0 - cnt/ny) + np.log(ymax)*(cnt/ny)) for cnt in range(ny)]

                if depth is None:
                    bins = np.zeros((ny-1,nx-1))
                else:
                    bins = np.zeros((depth,ny-1,nx-1))
                
                return bins, xgrid, ygrid

#def coord_bin_cent(i_m, i_r):
#   return [mgrid[i_m]+dy/2., rgrid[i_r]+dx/2.]

def t_lin(v,vi,vip1):
    return (v-vi)/(vip1-vi)

def binning_MR(R_data, M_data, rgrid, mgrid, stab, bins, mode=u.Mode.BOTH):

    dx = rgrid[1]-rgrid[0] 
    dy = mgrid[1]-mgrid[0]

    # Filter the data
    if mode == u.Mode.NS:
        branches_to_plot = [1] if np.max(stab) == 1 else []
    elif mode == u.Mode.TWINS:
        branches_to_plot = [1,2] if np.max(stab) == 2 else []
    else:
        branches_to_plot = list(range(1,int(np.max(stab))+1))

    for branch in branches_to_plot:

        R = R_data[np.where(stab==branch)]
        M = M_data[np.where(stab==branch)]

        if u.MMIN < np.max(M):
            ind = np.searchsorted(M,u.MMIN)

            i_m = np.searchsorted(mgrid, M[ind]) - 1
            i_r = np.searchsorted(rgrid, R[ind]) - 1
            
            bins[i_m,i_r] += 1

            for cnt in range(len(M[M>u.MMIN])-1):
                ind += 1

                dir_r = int(np.sign(R[ind] - R[ind-1]))
                dir_m = 1 #always growing

                #number of crossings
                n_m = int((M[ind]-mgrid[i_m])/dy)

                m_cross = [[t_lin(mgrid[ind_m], M[ind-1], M[ind]), ind_m, 'm', mgrid[ind_m]] for ind_m in (i_m + dir_m*(1+np.arange(n_m)))]
                if dir_r > 0:
                    n_r = int((R[ind] - rgrid[i_r])/dx) # int part
                    r_cross = [[t_lin(rgrid[ind_r], R[ind-1], R[ind]), ind_r, 'r', rgrid[ind_r]] for ind_r in (i_r + dir_r*(1+np.arange(n_r)))]
                else:
                    n_r = int((rgrid[i_r] + dx - R[ind])/dx) # int part
                    r_cross = [[t_lin(rgrid[ind_r], R[ind-1], R[ind]), ind_r-1, 'r', rgrid[ind_r]] for ind_r in (i_r + dir_r*(np.arange(n_r)))]
                crosses = sorted(m_cross+r_cross)
                    
                for t, index, m_or_r, val in crosses:
                    if m_or_r == 'r':
                        i_r = index
                        bins[i_m,i_r] += 1

                    if m_or_r == 'm':
                        i_m = index
                        bins[i_m,i_r] += 1
  
        else:
            pass 


def binning_EOS(x_data, y_data, xgrid, ygrid, stab, bins, mode=u.Mode.BOTH):
    """
    Bin EOS
    """
    
    # Filter
    if mode == u.Mode.NS and np.max(stab) != 1:
        return
    elif mode == u.Mode.TWINS and np.max(stab) !=2:
        return


    ind = np.searchsorted(x_data, u.EMINH)

    i_x = np.searchsorted(xgrid, x_data[ind]) - 1
    i_y = np.searchsorted(ygrid, y_data[ind]) - 1

    bins[i_y,i_x] += 1

    for cnt in range(np.searchsorted(x_data, u.EMAXH, side='left')-ind-1): 
        ind += 1

        dir_x = 1 #always growing
        dir_y = int(np.sign(y_data[ind] - y_data[ind-1]))
        

        #number of crossings
        n_x = np.searchsorted(xgrid, x_data[ind]) - np.searchsorted(xgrid, x_data[ind-1])
        n_y = np.searchsorted(ygrid, y_data[ind]) - np.searchsorted(ygrid, y_data[ind-1])

        x_cross = [[t_lin(xgrid[ind_x], x_data[ind-1], x_data[ind]), ind_x, 'x', xgrid[ind_x]] for ind_x in (i_x + dir_x*(1+np.arange(n_x)))]
        y_cross = [[t_lin(ygrid[ind_y], y_data[ind-1], y_data[ind]), ind_y, 'y', ygrid[ind_y]] for ind_y in (i_y + dir_y*(1+np.arange(n_y)))]

        crosses = sorted(y_cross+x_cross)

        for t, index, y_or_x, val in crosses:
            if y_or_x == 'x':
                i_x = index
                bins[i_y,i_x] += 1


            if y_or_x == 'y':
                i_y = index
                bins[i_y,i_x] += 1



def KDE(ax, xgrid, ygrid, bins):
    """
    Plot KDE
    """
    xcent = [(xgrid[i]+xgrid[i+1])/2 for i in range(len(xgrid)-1)]
    ycent = [(ygrid[i]+ygrid[i+1])/2 for i in range(len(ygrid)-1)]

    xmesh, ymesh = np.meshgrid(xcent, ycent)

    sns.kdeplot(data=pd.DataFrame(np.array(
        [xmesh.ravel(), ymesh.ravel(), bins.ravel()]).T, 
        columns=['x','y','w']),
        x='x', y='y', weights='w',
        cmap="Blues", 
        fill=False, 
        levels = [1-0.95,1-0.68], 
        bw_adjust = 0.25,
        ax = ax)



if __name__ == "__main__":

    import h5py
    import numpy as np
    import matplotlib.pyplot as plt
    from plot_settings import EOS_settings, subplot_settings
    from utils import nuclPoints, PRXboundary
    from classes import Style, Mode, parse_args
    #import argparse

    dir = "../build/"
    name = "Hile-run4-0.h5"
    NEOS = 50

    ptot_max = 81.56746479870479
    nrow_fig = 5
    ncol_fig = 1

    fig, axs = subplot_settings(nrow_fig, ncol_fig, figsize=(5, 4), wspace=0.1, hspace=0.1, sharex=True, sharey=True, squeeze=False)
    args = parse_args()
    mode = Mode(args.mode)
    style = Style(args.style)

    #for j in range(ncol_fig):
    #    for i in range(nrow_fig):
    #        EOS_settings(axs[i,j])

    with h5py.File(f"{dir}{name}", 'r') as file:
        for i in range(NEOS):
            tov_path = f"/{i}/TOV"
            params_path = f"/{i}/params"
            eos_path = f"/{i}/EOS"

            e_data = file[f"{eos_path}/e"][:]
            p_data = file[f"{eos_path}/p"][:]
            M_data = file[f"{tov_path}/M"][:]  
            R_data = file[f"{tov_path}/R"][:]
            stab = file[f"{tov_path}/stab"][:]
            ptot = file[f"{params_path}/ptot"][:][0]
            const_2M = file[f"{params_path}/pM"][:][0]
            pCET = file[f"{params_path}/pCET"][:][0]
            pGW = file[f"{params_path}/pGW"][:][0]
            dn = file[f"{params_path}/dn"][:][0]
            npt = file[f"{params_path}/nPTl"][:][0]
            e_cent = file[f"{tov_path}/e_cent"][:]


            if np.max(stab) == 2:
                ptot = 80#pGW*const_2M*pCET


                #M_1, M_2, first_branch, second_branch = find_max(M_data, stab, mode=mode.BOTH)
                #print(M_1, M_2)
                #M_1, first_branch = find_max(M_data, stab, mode=mode.NS)
                #print(M_1)
                #test plots with constraint evolution
                #plot_MR(axs, 0, R_data, M_data, stab, likelihood = ptot, tot_likelihood=ptot_max, style=style.POSTERIOR, mode=mode.BOTH)
                plot_MR_categories(axs, R_data, M_data, stab, likelihood = ptot, tot_likelihood=ptot_max, style=style, mode=mode, col = 0)

                #sort_and_plot_twins_prior_MR(R_data, M_data, M_1, M_2, first_branch, second_branch, axs, col=0)
                #plot_twins_posterior_EOS(e_data, p_data, M_1, M_2, first_branch, second_branch, ptot, ptot_max, axs, col=0)
                #plot_EOS(axs, 0, e_data, p_data, M_1, first_branch, M_2=M_2, second_branch=second_branch, likelihood = ptot, tot_likelihood=ptot_max, style=style.POSTERIOR, mode=mode.BOTH)
                
                #sort_and_scatter_twins_posterior(npt/0.16, dn, M_1, M_2, ptot[0], ptot_max, axs, col=0)
                #sort_and_scatter_twins_prior(npt/0.16, dn, M_1, M_2, axs, col=0)

                #print(M_data[stab==1], e_cent[first_branch])
                #print('edata', e_data[108:120])
                #print('pdata'   , p_data[108:120])

                #test plots with twin stability
                #find last point with interpolation and plot 

                #first_interp, second_interp = interpolate(e_cent, first_branch, e_data, p_data), interpolate(e_cent, second_branch, e_data, p_data)
                #plot_twins_stability(e_cent, e_data, p_data, M_1, M_2, first_branch, second_branch, 1, ptot, axs[0,0])
                
                #interp2 = np.interp(e_cent[second_branch], e_data, p_data)
                #print(len(interp), len(first_branch))
                #axs[0,0].plot(first_interp[0], first_interp[1], 'o', c='C'+str(i), markeredgecolor='black')
                #axs[0,0].plot(second_interp[0], second_interp[1], 'o', c='C'+str(i))
                

                #axs[0,0].fill(nuclPoints[:,0], nuclPoints[:,1], facecolor=settings['limit_color'], edgecolor=settings['limit_color'], zorder=4)
                #axs[0,0].plot(PRXboundary[:,0],PRXboundary[:,1], color='black', zorder=3)
                #print(M_1)
                #print(second_branch)
    
    plt.show()