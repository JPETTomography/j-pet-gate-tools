#!/usr/bin/python

import sys, getopt, os
import jinja2
__author__ ="pkopka"

#Script render Geometry macro for Gate simulation
#Create Gemerty.mac with <n_strip> strips  splited on <n_part> regions.
#Template pattern: Geometry_384_strips_L50.mac made by Pawel Kowalski.
#python splitter.py -l <strip_lenght> -n <part_number> -o <outputfile> -s <number_strips>

def main(argv):
    """
    Get options and argument from command line.
    :param argv:
    :return:
    """
    outputfile = "Geometry_strips_L50.mac"
    n_part = 50 # number of detector parts
    l = 50  # cm
    n_strip = 1 # nymber of strips
    try:
        opts, args = getopt.getopt(argv, "hl:n:o:s:", ["number-part=", "ofile=","strip-lenght=","strip-number="])
    except getopt.GetoptError:
        print('splitter.py -l <strip_lenght> -n <part_number> -o <outputfile> -s <number_strips>')
        sys.exit(2)
    for opt, arg in opts:
        if opt == '-h':
            print('splitter.py -l <strip_lenght> -n <part_number> -o <outputfile>')
            sys.exit()
        elif opt in ("-n"):
            try:
                n_part = int(arg)
            except ValueError:
                print("number of part must be integer")
        elif opt in ("-s"):
            try:
                n_strip = int(arg)
            except ValueError:
                print("number ofstrips must be integer")
        elif opt in ("-o"):
            outputfile = arg
        elif opt in ("-l"):
            try:
                l = float(arg)
            except ValueError:
                print("lenght of strip  must be float")
    print(n_part, l, n_strip, outputfile)
    return (n_part, l, n_strip, outputfile)


def render(tpl_path, context):
    """
    Render Geometry template.

    :param tpl_path:
    :param context:
    """
    path, filename = os.path.split(tpl_path)
    return jinja2.Environment(
        loader=jinja2.FileSystemLoader(path or './')
    ).get_template(filename).render(context)


if __name__ == "__main__":
    (n_part, l, n_strip, outputfile) = main(sys.argv[1:])
    context = dict(L=l, n_part=n_part, n_strip=n_strip, color=['blue', 'green'])
    path_file = os.path.dirname(os.path.realpath(__file__))
    result = render(os.path.join(path_file,'templates/Geometry_strips_L50.temp'), context)
    model_file = open(outputfile, 'w')
    model_file.write(result)
