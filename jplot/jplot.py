from mpl_toolkits.mplot3d import Axes3D
import matplotlib.pyplot as plot
import sys,os


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

if __name__ == "__main__":
    inputfile = main(sys.argv[1:])
    jplot(inputfile)