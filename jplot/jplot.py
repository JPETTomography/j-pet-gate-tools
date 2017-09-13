from mpl_toolkits.mplot3d import Axes3D
import matplotlib.pyplot as plot
import sys,os
import ROOT
import numpy as np
from matplotlib.colors import LogNorm

__author__ ="pkopka"

#Script plot line of response
# input file generate aplication analusis made by Pawel Kowalski
# python jplot.py <inputfile>

def main(argv):
    """
    Get argument from command line.
    :param argv:
    :return:
    """
    inputfile = "coincidences.txt"
    if argv:
        inputfile=argv[0]

    return inputfile




def jplot(inputfile):
    """
    Plot line of response from <inputfile> file with coincidences

    :param inputfile: X1 Y1 Z1 time1 X2 Y2 Z2 time2 ID1 ID2 Energy1 Energy2 Type
    :return:
    """

    fig = plot.figure()
    ax = Axes3D(fig)
    ax.set_xlabel('X axis')
    ax.set_ylabel('Y axis')
    ax.set_zlabel('Z axis')

    path_file = os.path.dirname(os.path.realpath(__file__))
    try:
        f = open(os.path.join(path_file,inputfile),'r')
    except IOError:
        print("No such file as %s" % os.path.join(path_file,inputfile))
        sys.exit(2)
    for line in f.readlines():
        line = list(map(float,line.split()))
        # print(line)
        try:
            X= [line[0],line[4]]
            Y = [line[1], line[5]]
            Z = [line[2], line[6]]
            if line[12] ==1.:
                color = 'red'
            else:
                color = 'blue'
        except IndexError:
            print("Bad structure. X1 Y1 Z1 time1 X2 Y2 Z2 time2 ID1 ID2 Energy1 Energy2 Type")
            sys.exit(2)
        ax.plot_wireframe(X, Y, Z, color=color)


    plot.show()

def jplot_TOF(inputfile):
    """
    Plot line of response from <inputfile> file with coincidences

    :param inputfile: X1 Y1 Z1 time1 X2 Y2 Z2 time2 ID1 ID2 Energy1 Energy2 Type
    :return:
    """

    fig = plot.figure()
    ax = Axes3D(fig)
    ax.set_xlabel('X axis')
    ax.set_ylabel('Y axis')
    ax.set_zlabel('Z axis')

    path_file = os.path.dirname(os.path.realpath(__file__))
    try:
        f = open(os.path.join(path_file,inputfile),'r')
    except IOError:
        print("No such file as %s" % os.path.join(path_file,inputfile))
        sys.exit(2)
    for line in f.readlines():
        line = list(map(float,line.split()))
        # print(line)
        try:
            X= [line[0],line[4]]
            Y = [line[1], line[5]]
            Z = [line[2], line[6]]
            if line[12] ==1.:
                color = 'red'
            else:
                color = 'blue'
        except IndexError:
            print("Bad structure. X1 Y1 Z1 time1 X2 Y2 Z2 time2 ID1 ID2 Energy1 Energy2 Type")
            sys.exit(2)
        ax.plot_wireframe(X, Y, Z, color=color)


    plot.show()

def jplot_diff(inputfile_1, inputfile_2 ):
    """
    Plot line of response from two lm inputfile on the same figure

    :param inputfile: X1 Y1 Z1 time1 X2 Y2 Z2 time2 ID1 ID2 Energy1 Energy2 Type
    :return:
    """

    fig = plot.figure()
    ax = Axes3D(fig)
    ax.set_xlabel('X axis')
    ax.set_ylabel('Y axis')
    ax.set_zlabel('Z axis')

    path_file = os.path.dirname(os.path.realpath(__file__))
    try:
        f = open(os.path.join(path_file,inputfile_1),'r')
    except IOError:
        print("No such file as %s" % os.path.join(path_file,inputfile_1))
        sys.exit(2)
    for line in f.readlines():
        line = list(map(float,line.split()))
        # print(line)
        try:
            X= [line[0],line[4]]
            Y = [line[1], line[5]]
            Z = [line[2], line[6]]

        except IndexError:
            print("Bad structure. X1 Y1 Z1 time1 X2 Y2 Z2 time2 ID1 ID2 Energy1 Energy2 Type")
            sys.exit(2)
        ax.plot_wireframe(X, Y, Z, color='green')


    try:
        f = open(os.path.join(path_file,inputfile_2),'r')
    except IOError:
        print("No such file as %s" % os.path.join(path_file,inputfile_2))
        sys.exit(2)
    for line in f.readlines():
        line = list(map(float,line.split()))
        # print(line)
        try:
            X= [line[0],line[4]]
            Y = [line[1], line[5]]
            Z = [line[2], line[6]]
        except IndexError:
            print("Bad structure. X1 Y1 Z1 time1 X2 Y2 Z2 time2 ID1 ID2 Energy1 Energy2 Type")
            sys.exit(2)
        ax.plot_wireframe(X, Y, Z, color='red')


    plot.show()




def jplot_cylindrical(inputfile_root,branch="Coincidences",r = 450,
               n_strip = 384.0,len_crystal=10, down =-240):
    """
    Plot LOR from root file. True Coincidences position and crystal position on the same figure.

    :param inputfile_root:
    :param branch:
    :param r: iternal radius of cylindrical scanner
    :param n_strip: number of strips
    :param len_crystal: Z size (mm)
    :param down: Z offset
    :return:
    """

    from math import cos, sin, pi
    f = ROOT.TFile(inputfile_root)
    print f.ls()
    t = f.Get(branch)
    fig = plot.figure()
    ax = Axes3D(fig)
    ax.set_xlabel('X axis')
    ax.set_ylabel('Y axis')
    ax.set_zlabel('Z axis')

    alfa = 2.0*pi/n_strip
    for event in t:
        X = [event.globalPosX1, event.globalPosX2]
        Y = [event.globalPosY1, event.globalPosY2]
        Z = [event.globalPosZ1, event.globalPosZ2]
        ax.plot_wireframe(X, Y, Z)

        X_2 = [r*cos(event.rsectorID1*alfa+alfa/2.0), r*cos(event.rsectorID2*alfa+alfa/2.0)]
        Y_2 = [r*sin(event.rsectorID1*alfa+alfa/2.0), r*sin(event.rsectorID2*alfa+alfa/2.0)]
        Z_2 = [event.crystalID1*len_crystal-len_crystal/2.0+down, event.crystalID2*len_crystal-len_crystal/2.0+down]
        ax.plot_wireframe(X_2, Y_2, Z_2, color='red')
    plot.show()

def spectrum(input_file,strips =384, t_size=192, s_size =192, t_window =3000, title =''):
    """
    Plot difference between rsectorIDs and times
    :param input_file:
    :return:
    """
    fig = plot.figure()
    path_file = os.path.dirname(os.path.realpath(__file__))
    try:
        f = open(os.path.join(path_file, input_file), 'r')
    except IOError:
        print("No such file as %s" % os.path.join(path_file, input_file))
        sys.exit(2)
    matrix = np.zeros((s_size,t_size))
    matrix.astype(int)
    diff_v = np.linspace(0, int(strips/2), num=s_size, endpoint=True)
    time_v= np.linspace(0,t_window,num=t_size,endpoint=True)
    time_list=[]
    diff =[]
    p=0
    l_num =0
    for line in f.readlines():
        line = list(map(float, line.split()))
        l_num+=1
        # try:
        if abs(line[3]-line[7]) >= t_window:

            p+=1
        else:
            try:

                diffid= np.where(int(min(abs(line[8]-line[9]),strips-abs(line[8]-line[9])))<=diff_v)
                time= np.where(abs(line[3]-line[7])<time_v)
                if not diffid[0].any():
                    i = s_size - 1
                else:
                    i = diffid[0][0]

                if not time[0].any():
                    j=t_size-1
                else:
                    j= time[0][0]
                # print 'diff ID',i
                # print 'time', j
                matrix[i,j]+=1
                # time_list.append(abs(line[3] - line[7]))
                # diff.append(min(abs(line[8] - line[9]), strips - abs(line[8] - line[9])))
            except KeyError:
            # except IndexError:
                matrix[i, t_size-1] += 1


        # except IndexError:
        #     print("Bad structure. X1 Y1 Z1 time1 X2 Y2 Z2 time2 ID1 ID2 Energy1 Energy2 Type")
        #     sys.exit(2)
    print 'p :', p
    print l_num
    print 'sum', sum(matrix.flatten()),'max:', max(matrix.flatten())

    imgplot = plot.imshow(np.fliplr(np.rot90(matrix,2)),extent=[0,300,0,int(strips/2)],norm=LogNorm(vmin=0.01,vmax=max(matrix.flatten())) )
    # imgplot.set_cmap('binary')
    plot.xlabel(r"$\Delta t [10^{-11} s]$")
    plot.ylabel(r"$\Delta$ rsectorID")
    plot.title(title)
    plot.colorbar()
    plot.savefig('../assets/%s.png' % title)
    # plot.figure()
    # plot.plot(time_list,diff,'o')


if __name__ == "__main__":
    inputfile = main(sys.argv[1:])
    jplot(inputfile)
    # jplot('../assets/lm_files/25sourceHE.txt')
    # for file_name in ['000', '0020', '2500', '4000']:
    #     for b in [ 'HE']:
    #         spectrum('../assets/lm_files/lm_cut_%s%s.txt' % (file_name,b),t_size=192,s_size=192,title='pozycja %s %s - prawdziwe' % (file_name, b[:2]))
    # spectrum('../assets/lm_files/con_25.txt',t_size=192,s_size=192)
    # spectrum('../assets/lm_files/con_25_50.txt',t_size=192,s_size=192)
    #     plot.show()
