"""
====================================================================
Linear Discriminant Analysis with covariance ellipsoid
====================================================================

This script plots the covariance ellipsoids of each class and
decision boundary learned by LDA. The ellipsoids display
the double standard deviation for each class. With LDA, the
standard deviation is the same for all the classes.

"""

# %%
# Colormap
# --------

import matplotlib.pyplot as plt
import matplotlib as mpl
import numpy as np                   ##for creating visualization
from matplotlib import colors              ## to specifie colors

cmap = colors.LinearSegmentedColormap(
    "red_blue_classes",
    {
        "red": [(0, 1, 1), (1, 0.7, 0.7)],
        "green": [(0, 0.7, 0.7), (1, 0.7, 0.7)],
        "blue": [(0, 0.7, 0.7), (1, 1, 1)],
    },
)
plt.cm.register_cmap(cmap=cmap)


def dataset_fixed_cov():
    """Generate 2 Gaussians samples with the same covariance matrix"""
    n, dim = 100, 2
    [Muons, Slope] = np.loadtxt("proton_hits_beta.txt", unpack = True)
    [Muons, Slope] = np.loadtxt("Iron_hits_beta.txt", unpack = True)

    X = np .r_[np.loadtxt("proton_hits_beta.txt"), np.loadtxt("Iron_hits_beta.txt")]

    y = np.hstack((np.zeros(n), np.ones(n)))
    return X, y

# Plot functions
# --------------

from scipy import linalg


def plot_data(lda, X, y, y_pred, fig_index):
    splot = plt.subplot(1, 1, fig_index)

    if fig_index == 1:
        plt.title("Proton vs Iron", fontsize=25)
        plt.xlabel("$dE_{Reco}$ (GeV)", fontsize=20) 

        plt.ylabel("\u03B2", fontsize=20)

    tp = y == y_pred  # True Positive
    tp0, tp1 = tp[y == 0], tp[y == 1]
    X0, X1 = X[y == 0], X[y == 1]
    X0_tp, X0_fp = X0[tp0], X0[~tp0]
    X1_tp, X1_fp = X1[tp1], X1[~tp1]

    # class 0: dots
    plt.scatter(X0_tp[:, 0], X0_tp[:, 1], marker=".", color="red", label="True H")
    plt.legend()
    plt.scatter(X0_fp[:, 0], X0_fp[:, 1], marker="x", s=20, color="#990000", label="False H")  # dark red
    plt.legend()
    # class 1: dots
    plt.scatter(X1_tp[:, 0], X1_tp[:, 1], marker=".", color="blue", label="True Fe")
    plt.legend()
    plt.scatter(
        X1_fp[:, 0], X1_fp[:, 1], marker="x", s=20, color="#000099", label="False Fe"
    )  # dark blue
    plt.legend()
    # class 0 and 1 : areas
    nx, ny = 200, 200
    x_min, x_max = plt.xlim()
    y_min, y_max = plt.ylim()
    xx, yy = np.meshgrid(np.linspace(x_min, x_max, nx), np.linspace(y_min, y_max, ny))
    Z = lda.predict_proba(np.c_[xx.ravel(), yy.ravel()])
    Z = Z[:, 1].reshape(xx.shape)
    plt.pcolormesh(
        xx, yy, Z, cmap="red_blue_classes", norm=colors.Normalize(0.0, 1.0), zorder=0
    )
    plt.contour(xx, yy, Z, [0.5], linewidths=2.0, colors="white")

    # means
    plt.plot(
        lda.means_[0][0],
        lda.means_[0][1],
        "*",
        color="yellow",
        markersize=15,
        markeredgecolor="grey",
    )
    plt.plot(
        lda.means_[1][0],
        lda.means_[1][1],
        "*",
        color="yellow",
        markersize=15,
        markeredgecolor="grey",
    )

    return splot


def plot_ellipse(splot, mean, cov, color):
    v, w = linalg.eigh(cov)
    u = w[0] / linalg.norm(w[0])
    angle = np.arctan(u[1] / u[0])
    angle = 180 * angle / np.pi  # convert to degrees
    # filled Gaussian at 2 standard deviation
    ell = mpl.patches.Ellipse(
        mean,
        2 * v[0] ** 0.5,
        2 * v[1] ** 0.5,
       180 + angle,
      facecolor=color,
      edgecolor="black",
        linewidth=2,
    )
    ell.set_clip_box(splot.bbox)
    ell.set_alpha(0.2)
    splot.add_artist(ell)
    splot.set_xticks(())
    splot.set_yticks(())


def plot_lda_cov(lda, splot):
    plot_ellipse(splot, lda.means_[0], lda.covariance_, "red")
    plot_ellipse(splot, lda.means_[1], lda.covariance_, "blue")


plt.figure(figsize=(10, 8), facecolor="white")
plt.suptitle(
    "",
    y=0.98,
    fontsize=20,
)

from sklearn.discriminant_analysis import LinearDiscriminantAnalysis


for i, (X, y) in enumerate([dataset_fixed_cov()]):
    # Linear Discriminant Analysis
    lda = LinearDiscriminantAnalysis(solver="svd", store_covariance=True)
    y_pred = lda.fit(X, y).predict(X)
    lda.fit_transform(X,y)
    lda.transform(X)
    X_testProton = lda.transform(X[y == 0])
    X_testIron = lda.transform(X[y == 1])
    print(X_testProton)
    print(X_testIron)
    splot = plot_data(lda, X, y, y_pred, fig_index=2 * i + 1)
    plot_lda_cov(lda, splot)
    plt.axis("tight")



plt.tight_layout()
plt.subplots_adjust(top=0.92)
plt.xticks(np.arange(0, 2000000, 500000))
##plt.xscale("log")
plt.xticks(fontsize=20)
plt.yticks(np.arange(2.27, 3.4, 0.2))
plt.yticks(fontsize=20)
plt.show()
