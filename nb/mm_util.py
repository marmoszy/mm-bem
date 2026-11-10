# utility functions for python mm-bem codes
# plot_msh, plot_lbs, print_lbs, tsv_read, mfrustum, osphere
# MM 9.11.2026

# --- plot mesh ---
def plot_mesh(v, e, fname, id, path):
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
    x_range = x_max - x_min
    y_range = y_max - y_min
    z_range = z_max - z_min
    max_range = max(x_range, y_range, z_range)
    x_mid = (x_max + x_min) / 2
    y_mid = (y_max + y_min) / 2
    z_mid = (z_max + z_min) / 2
    ax.set_xlim(x_mid - max_range / 2, x_mid + max_range / 2)
    ax.set_ylim(y_mid - max_range / 2, y_mid + max_range / 2)
    ax.set_zlim(z_mid - max_range / 2, z_mid + max_range / 2)
    ax.set_box_aspect((1, 1, 1))
    ax.view_init(elev=10, azim=20, vertical_axis='y')
    ax.set_xlabel('X')
    ax.set_ylabel('Y')
    ax.set_zlabel('Z')
    ax.set_title(r"%s, %dx%d"%(fname,len(v),len(e)))
    from matplotlib.ticker import MaxNLocator
    ax.xaxis.set_major_locator(MaxNLocator(nbins=4))
    ax.yaxis.set_major_locator(MaxNLocator(nbins=4))
    ax.zaxis.set_major_locator(MaxNLocator(nbins=4))
    plt.savefig(path+fname[:-4]+"_"+id+".svg", format="svg", bbox_inches="tight")
    plt.show()

# --- plot lbs ---
def plot_lbs(psc, lbs, fname, id, path):
    import matplotlib.pyplot as plt
    from numpy import pi, log10, linspace
    vmin, vmax = -60, -20
    fig, ax = plt.subplots(subplot_kw={'projection': 'polar'})
    colors = ["red","green","blue","orange","cyan","magenta"]
    if len(psc)<10:
        angles = linspace(0, 2 * pi, 360, endpoint=False)
        for i in range(len(psc)):
            values_db = 20 * log10([abs(p) for p in psc[i]])
            ax.plot(angles, values_db, color=colors[i%len(colors)], label=r"%gkHz/%g$\degree$/%.2fdB"%(lbs[i][0]/1e3,lbs[i][1],20*log10(abs(lbs[i][2]))))
            ax.scatter(lbs[i][1]*pi/180, 20*log10(abs(lbs[i][2])), s=50, edgecolor='none', facecolor=colors[i%len(colors)], alpha=1, zorder=3)
    else:
        grouped = {}
        for l in lbs:
            if l[0] not in grouped: grouped[l[0]] = []  
            grouped[l[0]].append(l)
        for i, (f0, lbs1) in enumerate(grouped.items()):
            angles = [r[1]*pi/180 for r in lbs1]
            values_db = 20 * log10([abs(r[2]) for r in lbs1])
            ax.plot(angles, values_db, color=colors[1+i], label=r"%gkHz"%(lbs1[-1][0]/1e3))
        ax.set_title(r"%s TS($\theta$)"%(fname), va='bottom')
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
    plt.savefig(path+fname[:-4]+"_"+id+".svg", format="svg", bbox_inches="tight")
    plt.show()

# --- print lbs ---
def print_lbs(psc, lbs, fname, n, m, id, path, tm=None):
    s = "# mm-bem-%s %s: %s (%dx%d)"%(id.split("_")[0],tm!=None and "(%.2fs)"%tm or "",fname,n,m)
    for j in range(len(psc)):
        f0, th, lbs0, c0 = lbs[j][0], lbs[j][1], lbs[j][2], lbs[j][3]
        if len(psc)<10:
            s += "\n# th0=%g ph0=0 f0=%g c0=%g lbs=%g (%g,%g)\n"%(th,f0,c0,abs(lbs0),lbs0.real,lbs0.imag)
            s += "\n".join(["%d\t%g"%(i,abs(psc[j][i])) for i in range(len(psc[j]))])
        else:
            s += "\n%g %g %g %g %g"%(th,f0,lbs0.real,lbs0.imag,abs(lbs0))
    print(s)
    with open(path+fname[:-4]+"_"+id+".txt", "w") as fp: print(s,file=fp)

# --- tsv reader ---
def tsv_read(fname):
    with open(fname) as fp: lines = fp.readlines()
    return [[float(n) for n in l.split()] for l in lines if not l.strip().startswith('#')]

# --- outline conversion ---
def every(a, n):
    return [e for i, e in enumerate(a) if (i + n - 1) % n == n - 1]
def mfrustum(xyzwh, N = 24, W = 1):
    from math import sin, cos, pi
    xyzwh = [list(row) for row in zip(*xyzwh)] # transpose
    vx, fx = [], []
    n, nW = 0, W or 1
    w = every(xyzwh[3], nW)
    h = every(xyzwh[4], nW)
    x = every(xyzwh[0], nW)
    y = every(xyzwh[1], nW)
    z = every(xyzwh[2], nW)
    e1 = every(xyzwh[5], nW) if len(xyzwh[5]) > 0 else [1] * len(w)
    e2 = every(xyzwh[6], nW) if len(xyzwh[6]) > 0 else [1] * len(w)
    for i in range(len(xyzwh[3])):
        w[i], h[i] = w[i]/2, h[i]/2
    for j in range(len(w)):
        if w[j] == 0: w[j] = 0.0001
        if h[j] == 0: h[j] = 0.0001
        if j == 0: continue
        mode = 0
        if j == 1: mode |= 1
        if j == len(w) - 1: mode |= 2
        for i in range(N):
            sa, ca = math.sin(2*pi*i/N), math.cos(2*pi*i/N)
            for k in range(j-1, j+1):  # k takes values: j-1, j
                ss = math.copysign(1, sa) * (abs(sa) ** e1[k])
                sc = math.copysign(1, ca) * (abs(ca) ** e1[k])
                a, b, c = 0, w[k]*ss, h[k]*sc
                vx.append([a+x[k], b+y[k], c+z[k]])
        vx.append([x[j-1], y[j-1], z[j-1]])
        vx.append([x[j], y[j], z[j]])
        for i in range(N):
            fx.append([n+2*i, n+2*i+1, n+(2*i+2)%(2*N)])
            fx.append([n+(2*i+2)%(2*N), n+2*i+1, n+(2*i+3)%(2*N)])
            if mode & 1:
                fx.append([n+2*i, n+(2*i+2)%(2*N), n+2*N])
            if mode & 2:
                fx.append([n+2*i+1, n+2*N+1, n+(2*i+3)%(2*N)])                
        n += 2*N+2
    return (vx, fx)

# --- outline sphere/ellipsoid generator ---
"""
def osphere(a, b=None, c=None, n=7):
    from math import sin, cos, pi
    if b == None: b = a
    if c == None: c = b
    return [[-a*cos(i*pi/(n-1)),0,0,2*b*sin(i*pi/(n-1)),2*c*sin(i*pi/(n-1))] for i in range(n)]
"""
import math
def get_dynamic_param(param, i, total_steps, default_value=1):
    if not isinstance(param, (list, tuple)): return param
    if len(param) == 0: return default_value
    if len(param) == 1: return param[0]
    progress = 1.0 - (i / (total_steps - 1))
    exact_index = progress * (len(param) - 1)
    index_low = math.floor(exact_index)
    index_high = min(index_low + 1, len(param) - 1)
    weight = exact_index - index_low
    return param[index_low] * (1-weight) + param[index_high] * weight
def osphere(r, n=14, alpha=0, beta=0, e1=1, e2=None):
    if isinstance(r, (list, tuple)):
        a, b, c = r[0]
        b = r[1] if len(r) > 1 else a
        c = r[2] if len(r) > 2 else b
    else:
        a = b = c = r
    if n < 3: n = 3
    if e2 is None: e2 = e1
    x = []
    for i in range(n):
        rad_a = (get_dynamic_param(alpha, i, n, 0) * math.pi) / 180
        rad_b = (get_dynamic_param(beta, i, n, 0) * math.pi) / 180
        cur_e1 = get_dynamic_param(e1, i, n, 1)
        cur_e2 = get_dynamic_param(e2, i, n, 1)
        base_x = -a * math.cos(i * math.pi / (n - 1))
        row = [
            base_x * math.cos(rad_a) * math.cos(rad_b),
            base_x * math.sin(rad_b),
            base_x * math.sin(rad_a),
            2 * b * math.sin(i * math.pi / (n - 1)),
            2 * c * math.sin(i * math.pi / (n - 1)),
            cur_e1,
            cur_e2
        ]
        x.append(row)
    return x
    #s = "".join(["\t".join(map(str, row)) + "\n" for row in x])
    #return s

