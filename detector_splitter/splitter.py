#!/usr/bin/python

import sys, os
import jinja2
import argparse

__author__ ="pkopka"

#Script render Geometry macro for Gate simulation
#Create Geomerty.mac with <n_strip> strips  splited on <n_part> regions.
#JPET Geometry: max number strips - 384 Z size strip 500 mm
#python splitter.py -c <crystals_number> -o <outputfile> -s <number_strips>

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
    parser = argparse.ArgumentParser()
    parser.add_argument('--crystals','-c', type =int,
                        help ='number of crystals in strip', default=50)
    parser.add_argument('--strips','-s', type=int,
                        help='number of strips', default=348)
    parser.add_argument('-o','--ofile',help='output file',
                        default='../assets/jpet_geometry/jpet_geometry_L50.mac')
    args = parser.parse_args()
    args_dict = (vars(args))

    context = dict( number_of_crystal=args_dict['crystals'], number_of_strips=args_dict['strips'],
                    size=500.0/args_dict['crystals'])
    path_file = os.path.dirname(os.path.realpath(__file__))
    result = render(os.path.join(path_file,'templates/JPET_500mm_Clinder_sys.temp'), context)
    model_file = open(args_dict['ofile'], 'w')
    model_file.write(result)
