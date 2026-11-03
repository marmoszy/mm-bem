# --- plot mesh
def plot_mesh(v,e,fname,fname2):
    import matplotlib.pyplot as plt
    from mpl_toolkits.mplot3d.art3d import Poly3DCollection
    from numpy import array
    tri = [[array(v[i]) for i in t] for t in e]
    fig = plt.figure()
    ax = fig.add_subplot(111, projection='3d')
    mesh_collection = Poly3DCollection(tri, alpha=0.6)
    mesh_collection.set_facecolor('white')
    mesh_collection.set_edgecolor('darkblue')
    mesh_collection.set_linewidth(0.2)
    ax.add_collection3d(mesh_collection)
    v_arr = array(v)
    x_min, x_max = v_arr[:, 0].min(), v_arr[:, 0].max()
    y_min, y_max = v_arr[:, 1].min(), v_arr[:, 1].max()
    z_min, z_max = v_arr[:, 2].min(), v_arr[:, 2].max()
    ax.set_xlim(x_min, x_max)
    ax.set_ylim(y_min, y_max)
    ax.set_zlim(z_min, z_max)
    x_range = x_max - x_min
    y_range = y_max - y_min
    z_range = z_max - z_min
    ax.set_box_aspect((x_range, y_range, z_range))
    ax.view_init(elev=10, azim=20, vertical_axis='y')
    ax.set_xlabel('X')
    ax.set_ylabel('Y')
    ax.set_zlabel('Z')
    ax.set_title(r"%s, %dx%d"%(fname,len(v),len(e)))
    from matplotlib.ticker import MaxNLocator
    ax.xaxis.set_major_locator(MaxNLocator(nbins=4))
    ax.yaxis.set_major_locator(MaxNLocator(nbins=4))
    ax.zaxis.set_major_locator(MaxNLocator(nbins=4))
    plt.savefig(fname2, format="svg", bbox_inches="tight")
    plt.show()

# --- plot solution
def plot_lbs(psc,lbs,fname,fname2):
    import matplotlib.pyplot as plt
    from numpy import pi, log10, linspace
    angles = linspace(0, 2 * pi, 360, endpoint=False)
    vmin, vmax = -60, -20
    fig, ax = plt.subplots(subplot_kw={'projection': 'polar'})
    colors = ["red","green","blue","orange","cyan","magenta"]
    for i in range(len(psc)):
        values_db = 20 * log10([abs(p) for p in psc[i]])
        ax.plot(angles, values_db, color=colors[i%len(colors)], label=r"%gkHz/%g$\degree$/%.2fdB"%(lbs[i][0]/1e3,lbs[i][1],20*log10(abs(lbs[i][2]))))
        ax.scatter(lbs[i][1]*pi/180, 20*log10(abs(lbs[i][2])), s=50, edgecolor='none', facecolor=colors[i%len(colors)], alpha=1, zorder=3)
    ax.legend(loc='upper left', bbox_to_anchor=(1.02, 1))
    ax.set_ylim(vmin, vmax)
    ax.set_rticks(linspace(vmin, vmax, 5))
    ax.set_xticks(linspace(0, 2*pi-pi/6, 12))  # 30deg
    ax.set_rlabel_position(-20)
    ax.grid(True)
    if len(lbs)==1:
        ax.set_title(r"%s, $f_0$=%g kHz, $\theta_{BS}$=%g$\degree$"%(fname,lbs[-1][0]/1e3,lbs[-1][1]), va='bottom')
        ax.text(-0.065, -0.065, r"$l_{BS}=$%6f"%(abs(lbs[-1][2])),
            fontsize=11,color='black',transform=ax.transAxes)
        ax.text(0.775, -0.065, "TS=%.2f dB"%(20*log10(abs(lbs[-1][2]))),
            fontsize=11,color='black',transform=ax.transAxes)
    plt.savefig(fname2, format="svg", bbox_inches="tight")
    plt.show()

# --- print
def print_txt(fname,n,m,psc,lbs,id,fname2):
    s = "# mm-bem-%s: %s (%dx%d)"%(id,fname,n,m)
    for j in range(len(psc)):
        f0, th, lbs0, c0 = lbs[j][0], lbs[j][1], lbs[j][2], lbs[j][3]
        s += "\n# th0=%g ph0=0 f0=%g c0=%g lbs=%g (%g,%g)\n"%(th,f0,c0,abs(lbs0),lbs0.real,lbs0.imag)
        s += "\n".join(["%d\t%g"%(i,abs(psc[j][i])) for i in range(len(psc[j]))])
    print(s)
    with open(fname2, "w") as fp: print(s,file=fp)

# --- demo
"""
from mm_plots import plot_mesh, plot_lbs, print_txt
plot_mesh(v, e, fname, "yft_taichi_1.svg")
plot_lbs(psc, lbs, fname, "yft_taichi_2.svg")
print_txt(fname, n, m, psc, lbs, "taichi", "yft_taichi_3.txt")
"""
