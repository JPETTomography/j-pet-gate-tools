import sys
import LayersController

# =============================
if __name__ == '__main__':
    try:
        file_name = sys.argv[1]
    except IndexError:
        sys.exit("\nUsage: python reduce_multilayer.py <input file>  \n"
                 "\t [-o <output file>] [-l <layer ID>] \n"
                 "\t [-s <no of strips>]  [-with_tof] ")
    sys.stdout.write("-------------------------------------\n")
    sys.stdout.flush()

    # transformation options
    options = sys.argv[2:]

    # output filename
    if '-o' in options:
        out_filename = options[options.index('-o') + 1]
    else:
        out_filename = None

    # (OPTIONAL) all to lowercase
    options = map(lambda x: x.lower(), options)

    # check whether to correct TOF
    if '-with_tof' in options:
        lc = LayersController.LayersController(adjust_times=True)
    else:
        lc = LayersController.LayersController()

    # add 'ideal' radius 43.73 cm
    lc.add_layer(43.73)

    # default no of strips
    strips_in_total = 192
    # check layer if other than 1
    if '-l' in options:
        layer_id = int(options[options.index('-l') + 1])
        if layer_id == 0:
            # check no of strips
            if '-s' in options:
                strips_in_total = int(options[options.index('-s') + 1])
            lc.add_zero_layer(strips_in_total)
    else:
        layer_id = 1
    
    # print out the radius on a new single layer and the no of strips
    if layer_id == 0:
        radius_to_map = lc.geometry['ZeroLayerRadius']
    else:
	radius_to_map = lc.geometry['LayersRadiuses'][layer_id - 1]
    sys.stdout.write("The data will be mapped to radius R = %.2f cm " 
		     % radius_to_map)
    sys.stdout.write("(Number of strips: %s).\n\n" % strips_in_total)
    sys.stdout.flush()

    # import data
    lc.read_data(file_name)
    # !!! Here is the remap itself (static_radiuses could be set here) !!!
    lc.remap_data(out_layer_id=layer_id, static_radiuses=False)
    # export output
    lc.export_remapped_data(out_filename)


