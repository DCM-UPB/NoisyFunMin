from numpy import *
from matplotlib.pyplot import *
import matplotlib.animation as animation
import matplotlib.colors as colors

def getData(file):
    x = []
    y = []
    z = []
    for line in open(file):
        llist = line.split()
        if len(llist) == 6:  # could be valid x line
            if (llist[0] == 'x0' and llist[3] == 'x1'):  # yes it is
                x.append(float(llist[2]))
                y.append(float(llist[5]))
        if len(llist) == 5:  # could be valid f line
            if (llist[0] == 'f'):
                z.append(float(llist[2]))

    return [array(x), array(y), array(z)]


# Rosenbrock function
def rbf(x, y):
    return 100. * (y - x ** 2) ** 2 + (1. - x) ** 2


# showFigure script with hardcoded modes
ncg_split = 31


def showFigure(files, names, title, mode):  # mode can be 0 = cg, 1 = sgd, 2 = noisy-split, 3 = noisy-adam
    data = []
    for file in files:
        px, py, pz = getData(file)
        data.append([px, py, pz])

    # generate plot of Rosenbrock function
    x = linspace(-2., 2., 250)
    y = linspace(-1, 3, 250)
    X, Y = meshgrid(x, y)
    Z = rbf(X, Y)
    z_max = abs(Z).max()

    fig = figure()
    ax = fig.add_subplot(111)
    c = ax.pcolormesh(X, Y, Z, norm=colors.LogNorm(vmin=Z.min(), vmax=Z.max()), cmap='viridis')
    ax.set_title(title)
    ax.axis([x.min(), x.max(), y.min(), y.max()])
    fig.colorbar(c, ax=ax)

    clist = ['white', 'purple', 'pink', 'black', 'lime', 'red']
    mlist = ['o', 'x', 's', 'd', '*', '+']

    lines = []
    for it in range(len(files)):
        nlines = 1
        if mode == 2:  # we make split lines
            nlines = 2
        for il in range(nlines):
            line, = ax.plot([], [], '-', lw=1, color=clist[it + il])
            lines.append(line)
    ax.legend(names, loc='lower right')

    points = []
    for it in range(len(files)):
        npoints = 1
        if mode == 2:  # we make split lines
            npoints = 2
        for ip in range(npoints):
            point, = ax.plot([], [], marker=mlist[it + ip], color=clist[it + ip])
            points.append(point)

    time_template = 'step = %i'
    time_text = ax.text(0.05, 0.9, '', transform=ax.transAxes)

    def getReturnList(lines, points, time_text):  # assuming lines and points have equal len
        # I didn't find a better way of doing this
        if len(lines) == 0:
            return time_text
        if len(lines) == 1:
            return lines[0], points[0], time_text
        if len(lines) == 2:
            return lines[0], lines[1], points[0], points[1], time_text
        if len(lines) == 3:
            return lines[0], lines[1], lines[2], points[0], points[1], points[2], time_text
        if len(lines) == 4:
            return lines[0], lines[1], lines[2], lines[3], points[0], points[1], points[2], points[3], time_text
        if len(lines) == 5:
            return lines[0], lines[1], lines[2], lines[3], lines[4], points[0], points[1], points[2], points[3], points[
                4], time_text
        if len(lines) == 6:
            return lines[0], lines[1], lines[2], lines[3], lines[4], lines[5], points[0], points[1], points[2], points[
                3], points[4], points[5], time_text
        print("Only up to 6 lines supported.")

    def getTrueLen(n):
        if (mode != 2):
            return n
        else:
            if n > ncg_split:
                return ncg_split + (n - ncg_split) // 5
            else:
                return n

    def getTrueIdx(i):
        if (mode != 2):
            return i
        else:
            if i > ncg_split:
                return ncg_split + 5 * (i - ncg_split)
            else:
                return i

    def init():
        for line in lines:
            line.set_data([], [])
        for point in points:
            point.set_data([], [])
        time_text.set_text('')
        return getReturnList(lines, points, time_text)

    def animate(i):
        idx = getTrueIdx(i)
        if (mode != 2):  # split mode will get special treatment
            if (i > 0):
                for itl, line in enumerate(lines):
                    line.set_data(data[itl][0][0:idx + 1], data[itl][1][0:idx + 1])
            else:
                for itl, line in enumerate(lines):
                    line.set_data([data[itl][0][0], data[itl][0][0]], [data[itl][1][0], data[itl][1][0]])
        else:
            if (i > 0):
                for itl in range(len(lines) // 2):
                    if (i > ncg_split):
                        lines[itl + 1].set_data(data[itl][0][ncg_split:idx + 1], data[itl][1][ncg_split:idx + 1])
                    else:
                        lines[itl].set_data(data[itl][0][0:idx + 1], data[itl][1][0:idx + 1])
            else:
                for itl in range(len(lines) // 2):
                    lines[itl].set_data([data[itl][0][0], data[itl][0][0]], [data[itl][1][0], data[itl][1][0]])

        if (mode != 2):
            for itp, point in enumerate(points):
                point.set_data([data[itp][0][idx], data[itp][1][idx]])
        else:
            itshift = 0
            if (i > ncg_split):
                itshift = 1
            for itp in range(len(points) // 2):
                points[itp + itshift].set_data([data[itp][0][idx], data[itp][1][idx]])

        time_text.set_text(time_template % idx)

        return getReturnList(lines, points, time_text)

    def getInterval():
        if mode == 0:
            return 200
        if mode == 1:
            return 5
        if mode == 2:
            return 100
        if mode == 3:
            return 15

    ani = animation.FuncAnimation(fig, animate, arange(0, getTrueLen(len(px))), interval=getInterval(), blit=True,
                                  init_func=init)
    #    ani.save("anim_"+str(mode)+".mp4", fps=1000./getInterval())
    show()


# --- Script

prefix = "../../build/examples/"
showFigure([prefix + "cgsd.out", prefix + "cgfr.out", prefix + "cgpr.out", prefix + "cgpr0.out"],
           ["SD", "CG(FR)", "CG(PR)", "CG(PR0)"], "SD/CG variants, no noise", 0)
showFigure([prefix + "sgdm.out", prefix + "nest.out", prefix + "rmsp.out", prefix + "adad.out", prefix + "adam.out",
            prefix + "fire.out"], ["SGDM", "Nesterov", "RMSProp", "AdaDelta", "Adam", "FIRE"],
           "SGD algorithms (and FIRE), no noise", 1)
showFigure([prefix + "cg-sgd_noise.out"], ["CG", "SGDM"], "Noisy CG, followed by momentum SGD, with noise", 2)
showFigure(
    [prefix + "adam_noise.out", prefix + "fire_noise.out", prefix + "irene_noise.out", prefix + "irene2_noise.out"],
    ["Adam", "FIRE", "IRENE", "IRENE-B"], "Adam, FIRE and own variant IRENE), with noise", 3)
