import numpy as np
from mpl_toolkits.axes_grid1 import make_axes_locatable
import matplotlib.pyplot as plt
import h5py as h5

# set up figure
plt.rc("font", family="serif")
plt.rc("xtick", labelsize=7)
plt.rc("ytick", labelsize=7)
plt.rc("axes", labelsize=7)

#fig = plt.figure(figsize=(3.54, 6.0))  # single column fig
fig = plt.figure(figsize=(7.48, 6.0))  #two column figure
gs = plt.GridSpec(2, 2, wspace=0.0, hspace=0.35)

axs = []
axs.append(plt.subplot(gs[0, 0]))
axs.append(plt.subplot(gs[1, 0]))
axs.append(plt.subplot(gs[0, 1]))
axs.append(plt.subplot(gs[1, 1]))

for ax in axs:
    ax.minorticks_on()

axs[0].set_xlabel(r"$\Lambda_1$")
axs[0].set_ylabel(r"$\Lambda_1$")
axs[1].set_xlabel(r"$\Lambda_1$")
axs[1].set_ylabel(r"$\Lambda_1$")

axs[2].set_xlabel(r"$\Lambda$")
axs[3].set_xlabel(r"$\Lambda$")

axs[2].set_ylabel(r"pdf")
axs[3].set_ylabel(r"pdf")

axs[0].set_title(r"Original L/V data")
axs[1].set_title(r"Gaussian Mixture Model")


file_name = "low_spin_PhenomPNRT_posterior_samples.dat"
#file_name = "high_spin_PhenomPNRT_posterior_samples.dat"

# 0 costheta_jn
# 1 luminosity_distance_Mpc
# 2 m1_detector_frame_Msun
# 3 m2_detector_frame_Msun
# 4 lambda1
# 5 lambda2
# 6 spin1
# 7 spin2
# 8 costilt1
# 9 costilt2

data = np.loadtxt(file_name, skiprows=1)
print(np.shape(data))

xedges = np.linspace(0.0, 1600, 100)
yedges = np.linspace(0.0, 1600, 100)

#N = 10
#edges0 = np.linspace(np.min(data[:,0]), np.max(data[:,0]), N)
#edges1 = np.linspace(np.min(data[:,1]), np.max(data[:,1]), N)
#edges2 = np.linspace(np.min(data[:,2]), np.max(data[:,2]), N)
#edges3 = np.linspace(np.min(data[:,3]), np.max(data[:,3]), N)
#edges4 = np.linspace(np.min(data[:,4]), np.max(data[:,4]), N)
#edges5 = np.linspace(np.min(data[:,5]), np.max(data[:,5]), N)
#edges6 = np.linspace(np.min(data[:,6]), np.max(data[:,6]), N)
#edges7 = np.linspace(np.min(data[:,7]), np.max(data[:,7]), N)
#edges8 = np.linspace(np.min(data[:,8]), np.max(data[:,8]), N)
#edges9 = np.linspace(np.min(data[:,9]), np.max(data[:,9]), N)


H, xedges, yedges = np.histogram2d(data[:, 4], data[:, 5], bins=(xedges, yedges))

H = H.T  # Let each row list bins with common y range.

axs[0].imshow(
    H,
    interpolation="nearest",
    origin="low",
    extent=[xedges[0], xedges[-1], yedges[0], yedges[-1]],
)

H = H/np.sum(H)
axs[2].plot(xedges[:-1], np.diagonal(H))


# im = axs[0].imshow(
#            hdata_masked.T,
#            extent=[rad_grid[0], rad_grid[-1], mass_grid[0], mass_grid[-1] ],
#            origin='lower',
#            interpolation='nearest',
#            cmap="Reds",
#            vmin=0.0,
#            vmax=1.0,
#            aspect='auto',
#            )

# component number
if False:
    from sklearn.mixture import GaussianMixture as GMM

    data2d = np.vstack([data[:, 4], data[:, 5]]).T

    n_components = np.arange(1, 31)
    models = [
        GMM(n, covariance_type="full", random_state=0).fit(data2d) for n in n_components
    ]

    axs[1].plot(n_components, [m.bic(data2d) for m in models], label="BIC")
    axs[1].plot(n_components, [m.aic(data2d) for m in models], label="AIC")
    # axs[1].legend(loc='best')
    # axs[1].xlabel('n_components');


# 2D reduction
if True:
    # from sklearn.mixture import GMM
    from sklearn.mixture import GaussianMixture as GMM

    # reduce to L_1 / L_2 data only
    data2d = np.vstack([data[:, 4], data[:, 5]]).T

    gmm = GMM(n_components=10, covariance_type="full", random_state=0)
    gmm.fit(data2d)

    # display predicted scores by the model as a contour plot
    X, Y = np.meshgrid(xedges, yedges)
    XX = np.array([X.ravel(), Y.ravel()]).T
    Z = -gmm.score_samples(XX)
    Z = Z.reshape(X.shape)

    print(np.shape(Z))

    # CS = axs[1].contour(X, Y, Z,
    #        vmin=np.min(Z),
    #        vmax=np.max(Z),
    #        #norm=LogNorm(vmin=1.0, vmax=1000.0),
    #        #levels=np.logspace(0, 3, 10)
    #        )

    # transform back
    Z = np.exp(-Z)
    # Z = 10.0**(-1.0*Z)
    # Z = np.log(Z)
    print("min {} max {} ".format(np.min(Z), np.max(Z)))

    axs[1].imshow(
        Z,
        interpolation="nearest",
        origin="low",
        extent=[xedges[0], xedges[-1], yedges[0], yedges[-1]],
        cmap="viridis",
    )


    #get oneD
    Z[:,:] = Z[:,:]/np.sum(Z)
    oneD = np.diagonal(Z)

    #left, bottom, width, height = [0.55, 0.30, 0.2, 0.25]
    #ax2 = fig.add_axes([left, bottom, width, height])
    #ax2.plot(xedges, oneD)

    axs[3].plot(xedges, oneD)

    f5 = h5.File('LV_prior.h5', 'w')
    dset1 = f5.create_dataset("x",     data=xedges)
    dset2 = f5.create_dataset("y",     data=yedges)
    dset3 = f5.create_dataset("data",  data=Z)
    f5.close()


# Full 10D reconstruction
if False:
    from sklearn.mixture import GaussianMixture as GMM

    gmm = GMM(n_components=10, covariance_type="full", random_state=0)
    gmm.fit(data)

    # display predicted scores by the model as a contour plot
    X0, X1, X2, X3, X4, X5, X6, X7, X8, X9 = np.meshgrid(
            edges0,
            edges1,
            edges2,
            edges3,
            edges4,
            edges5,
            edges6,
            edges7,
            edges8,
            edges9)

    XX = np.array(
            [
        X0.ravel(), 
        X1.ravel(), 
        X2.ravel(), 
        X3.ravel(), 
        X4.ravel(), 
        X5.ravel(), 
        X6.ravel(), 
        X7.ravel(), 
        X8.ravel(), 
        X9.ravel(), 
           ]).T

    Z = -gmm.score_samples(XX)
    Z = Z.reshape(X.shape)

    print(np.shape(Z))


    # transform back
    Z = np.exp(-Z)
    print("min {} max {} ".format(np.min(Z), np.max(Z)))

    axs[1].imshow(
        Z,
        interpolation="nearest",
        origin="low",
        extent=[xedges[0], xedges[-1], yedges[0], yedges[-1]],
        cmap="viridis",
    )


# --------------------------------------------------
# GP estimate original full distribution
if False:
    from sklearn.gaussian_process import GaussianProcessRegressor
    from sklearn.gaussian_process.kernels import ConstantKernel, RBF

    # kernel
    rbf = ConstantKernel(1.0) * RBF(length_scale=1.0)
    #
    ##regressor
    gpr = GaussianProcessRegressor(
        kernel=rbf,
        #        alpha=noise**2,
    )

    ni, nf = np.shape(data)
    print("ni = {} nf = {}".format(ni, nf))
    target = np.ones((ni))

    gpr.fit(data, target)

    # mu_s, cov_s = gpr.predict(X, return_cov=True)


plt.savefig("LV_lambdas.pdf")
