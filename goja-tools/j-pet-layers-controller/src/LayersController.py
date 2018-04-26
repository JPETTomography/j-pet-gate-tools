import sys
import time
import numpy
import pandas


class LayersController:
    """
    Transforms 3-layer data into 1-layer data
    IMPORTANT: adjusting TOF (according to the change of layers) is set 
    as an instance parameter during __init__() as :param adjust_times, 
    while the mechanism to estimate new position, denoted in methods as :param static_radiuses,
    is 'hidden' and not set during initialisation (see the description in methods). 
    """
    # Class attribute
    C_LIGHT_SPEED = 299792458.  # m/s, float

    def __init__(self,
                 geometry=None,
                 adjust_times=False):
        """
        :param geometry: dictionary of parameters:
        :param adjust_times: whether to adjust times of hit or no  
        """
        if geometry is None:
            self.geometry = {
                'LayersRadiuses': [42.5, 46.75, 57.5],
                'ScannerLength': 50.,
                'StripCrossSection': [.7, 1.9],
                'ZeroLayerRadius': None
            }
        else:
            self.geometry = geometry

        # Internal variables (private)
        self.__lor_data = None
        self.__lor_output = None
        self.__adjust_times = adjust_times

        # Data state (Not private, since ends with '__')
        self.__data_state__ = {
            'filename': None,
            'loaded': False,
            'transformed': False
        }

    def add_zero_layer(self, total_strips_no=192):
        """
        Detects the radius which the strips are tightly composed in
        :param total_strips_no: number of strips (default is 48+48+96)
        :return: None
        """
        # (0.7/(2 sin(dQ/2))) * cos(dQ/2) + 0.95
        delta_theta = 2 * numpy.pi / total_strips_no
        self.geometry['ZeroLayerRadius'] = \
            round(self.geometry['StripCrossSection'][0] /
                  (2 * numpy.tan(delta_theta / 2)), 2) + self.geometry['StripCrossSection'][1] / 2

    def add_layer(self, radius):
        """
        Adds new radius to self.geometry['LayersRadiuses']
        :param radius: the radius of a new layer
        :return: None
        """
        self.geometry['LayersRadiuses'].append(radius)

    def read_data(self, filename):
        """
        Imports data from file in list mode
        :param filename: path to file
        :return: none
        """
        # state changed
        self.__data_state__['filename'] = filename
        sys.stdout.write("Loading data... ")
        self.__lor_data = pandas.read_table(filename, header=None)
        sys.stdout.write("Done!\nData has been successfully imported.\n")
        sys.stdout.flush()
        # state changed
        self.__data_state__['loaded'] = True

    def remap_data(self,
                   data=None,
                   out_layer_id=1,
                   static_radiuses=False,
                   adjust_times=None):
        """
        Transforms all imported data to single layer
        :param data: input data (pandas dataframe)
        :param out_layer_id: ID of the layer to coerce data
        :param adjust_times: whether to adjust times of hit or no
        :param static_radiuses: if False, the adjustment is made by comparison of hit positions,
                                otherwise the fixed fraction is used (R_new_layer/R_ini_layer)
        :return: message
        """
        # validate and set the unknown arguments
        if not self.__data_state__['loaded']:
            print "ERROR: No data imported yet. Please, import the data first."
            return "Aborted."
        if data is None:
            data = self.__lor_data
        if adjust_times is None:
            adjust_times = self.__adjust_times
        # create empty array
        self.__lor_output = numpy.zeros((self.__lor_data.shape[0], 8), dtype=float)
        # iterator
        iterator = 0
        # engage!
        sys.stdout.write('Mapping data:')  # setup progress bar
        for index, row in data.iterrows():
            self.__lor_output[iterator] = self.__event_to_layer(list(row),
                                                                out_layer_id,
                                                                static_radiuses,
                                                                adjust_times)
            iterator += 1
            step = round(float(data.shape[0]) / 100)  # 100 dots of progress
            if index % step == 0:
                sys.stdout.write('.')  # progress dots
                sys.stdout.flush()

        sys.stdout.write(" Done!\n")
        time.sleep(0.4)
        sys.stdout.flush()
        # state changed
        self.__data_state__['transformed'] = True
        return "Successfully finished the transformation"

    def export_remapped_data(self,
                             filename=None,
                             fmt='%.2f',
                             sep='\t'):
        """
        Exports __lor_output to a file
        :param filename: filename of output data in ASCII 
        :param fmt: output format (default is float rounded to 2 decimals)
        :param sep: delimiter for numpy.savetxt()
        :return: message
        """
        # validate if transformed
        if not self.__data_state__['transformed']:
            print "ERROR: No transformed data found."
            return "Aborted."
        # set the unset parameters
        if filename is None:
            filename = str(self.__data_state__['filename']) + '_remapped'
            # add a notation denoting TOF remapped
            if self.__adjust_times:
                filename += '_wTOF'
        # export data
        numpy.savetxt(filename, self.__lor_output, fmt=fmt, delimiter=sep)
        print "Successfully saved to file " + filename + "."
        return "Done!"

    # -------------- static methods ----------------
    @staticmethod
    def detect_lor_centre(lor):
        """
        Detect Cartesian coordinates for the centre of remapped LOR (for displacement calc)
        :param lor: line of response - single row from data
        :return: Cartesian triple
        """
        # 2D points (XY pair)
        p1_XY, p2_XY = lor[0:2], lor[4:6]
        if p1_XY[0] - p2_XY[0] == 0:
            XY = [p1_XY[0], 0.]
        elif p1_XY[1] - p2_XY[1] == 0:
            XY = [0., p1_XY[1]]
        else:
            # from line equation on a plane
            a = (p2_XY[1] - p1_XY[1]) / (p2_XY[0] - p1_XY[0])
            b = (p2_XY[0] * p1_XY[1] - p1_XY[0] * p2_XY[1]) / (p2_XY[0] - p1_XY[0])
            # the slope for the perpendicular is inv minus slope of the original (a -> -1/a)
            x = -(a * b) / (a ** 2 + 1)
            y = -x / a
            XY = [x, y]
        # fractions along original LOR from the 'centre' of remapped LOR
        norm_radiuses = [numpy.sqrt((lor[0] - XY[0]) ** 2 + (lor[1] - XY[1]) ** 2),
                         numpy.sqrt((lor[4] - XY[0]) ** 2 + (lor[5] - XY[1]) ** 2)]
        # from the similarity of triangles
        Z = (norm_radiuses[0] * lor[6] + norm_radiuses[1] * lor[2]) / sum(norm_radiuses)
        # output
        return [XY[0], XY[1], Z]

    @staticmethod
    def map_by_radius(cartesians,
                      r_new,
                      r_ini=None):
        """
        Remaps Cartesians using the relation between radiuses (similar triangles)   
        :param cartesians: list of [x, y, z]
        :param r_new: new radius
        :param r_ini: initial radius
        :return: list of [x, y, z]
        """
        if r_ini is None:
            r_ini = numpy.sqrt(cartesians[0] ** 2 + cartesians[1] ** 2)
        map_coefficient = r_new / r_ini
        return map(lambda x: x * map_coefficient, cartesians)

    # --------- private methods --------------------

    def __detect_layer(self, hit):
        """
        Estimates the layer at which the hit has happened 
        :param hit: one hit (at least X and Y)
        :return: layer radius, None if exceeds ranges
        """
        hit_radius = numpy.sqrt(hit[0] ** 2 + hit[1] ** 2)
        is_in_layer = map(lambda x:
                          x - self.geometry['StripCrossSection'][1] / 2
                          <= hit_radius
                          <= x + self.geometry['StripCrossSection'][1] / 2,
                          self.geometry['LayersRadiuses'])
        try:
            output = self.geometry['LayersRadiuses'][is_in_layer.index(True)]
        except ValueError:
            output = None
        return output

    def __hit_to_layer(self,
                       hit,
                       lor_centre,
                       out_layer_id=1,
                       static_radiuses=False,
                       adjust_times=False):
        """
        Remaps hit coordinates to new layer
        :param hit: coordinates + time tag
        :param lor_centre: where displacement (perpendicular) intercepts the LOR 
        :param out_layer_id: ID of the layer to coerce data
        :param static_radiuses: if False, the adjustment is made by comparison of hit positions,
                                otherwise the fixed fraction is used (R_new_layer/R_ini_layer)
        :param adjust_times: whether to adjust times of hit or no
        :return: new hit (coordinates + time tag)
        """
        lor_displacement_squarred = lor_centre[0] ** 2 + lor_centre[1] ** 2
        ini_layer = self.__detect_layer(hit)
        # outside strips - leave as it is
        if ini_layer is None:
            return hit
        # no action if the layers are the same
        if ini_layer == out_layer_id:
            return hit
        else:
            # not using numpy for vector subtraction is actually faster
            # biased hit corresponds to the move of (0,0,0) to lor_centre
            hit_biased = map(sum, zip(hit[0:3], map(lambda x: -x, lor_centre)))
            # hit distance from the centre of coordinates (in XY)
            hit_XY_squarred = hit[0] ** 2 + hit[1] ** 2
            # the difference between actual scattering ang geometrical centre of the strip
            depth_bias = numpy.sqrt(hit_XY_squarred) - ini_layer

            # choose the layer to map
            if out_layer_id == 0:
                # note that it has to be already defined (not None)!
                new_layer = self.geometry['ZeroLayerRadius']
            else:
                new_layer = self.geometry['LayersRadiuses'][out_layer_id - 1]

            # detect R_new, R_ini (for map_by_radius)
            if static_radiuses:
                r_new = numpy.sqrt(new_layer ** 2 - lor_displacement_squarred)
                r_ini = numpy.sqrt(ini_layer ** 2 - lor_displacement_squarred)
            else:
                r_new = numpy.sqrt((new_layer + depth_bias) ** 2 - lor_displacement_squarred)
                r_ini = numpy.sqrt(hit_XY_squarred - lor_displacement_squarred)

            # estimate new coordinates
            new_coordinates = self.map_by_radius(hit_biased, r_new, r_ini)
            # now move the centre of coordinates back to (0,0,0)
            output = map(sum, zip(new_coordinates, lor_centre))
            # adjust time if needed
            if adjust_times:
                # difference between LOR's propagation distance
                delta_distance = numpy.sqrt(sum(map(lambda x: x ** 2, new_coordinates))) \
                                 - numpy.sqrt(sum(map(lambda x: x ** 2, hit_biased)))
                # multiply by 100 to convert m to cm and by 1e12 to convert to ps
                # a warning should be ignored
                output.append(hit[3] + round(1e12 * delta_distance / (self.C_LIGHT_SPEED * 100.), 1))
            else:
                output.append(hit[3])
            return output

    def __event_to_layer(self,
                         event,
                         out_layer_id=1,
                         static_radiuses=False,
                         adjust_times=False):
        """
        Remaps hit coordinates (event=pair of hits) to new layer
        :param event: pair of hits (coordinates + time tags)
        :param out_layer_id: ID of the layer to coerce data
        :param adjust_times: whether to adjust times of hit or no
        :param static_radiuses: if False, the adjustment is made by comparison of hit positions,
                                otherwise the fixed fraction is used (R_new_layer/R_ini_layer)
        :return: new event (coordinates + time tags)
        """
        # preset parameters
        lor_centre = self.detect_lor_centre(event)
        # engage!
        remapped_event = map(lambda x: self.__hit_to_layer(x,
                                                           lor_centre,
                                                           out_layer_id,
                                                           static_radiuses,
                                                           adjust_times), [event[:4], event[4:8]])
        # flatten list
        return [item for sublist in remapped_event for item in sublist]
