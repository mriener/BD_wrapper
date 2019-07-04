import itertools
import os

import numpy as np

from astropy.io import fits
from astropy.table import Table


class BayesianDistance(object):
    def __init__(self, filename=None):
        """
        initializes the BayesianDistance class

        Parameters
        ----------
        path_to_bde : file path to the Bayesian distance program
        bde_script : read in fortran script of the Bayesian distance
            estimator
        path_to_file : file path to weighted FITS cube containing information
            about the decomposed Gaussians
        path_to_table : file path of the astropy table which contains the
            distance results
        path_to_input_table : file path of a table containing information of the
            Gaussian decompositions
        input_table : table containing information of the Gaussian
            decompositions
        save_input_table : The default is 'False'. If set to 'True' `input_table`
            is saved in the directory `path_to_table`
        verbose : The default is 'True'. Prints status messages to the
            terminal.
        gpy_setting : The default is 'False'. Set it to 'True' if `input_table`
            is created from a GaussPy decomposition
        intensity_threshold : Sets the threshold in integrated intensity of
            which decomposed Gaussian components should be considered. The
            default is '0.1'.
        distance_spacing : Only used for the creation of ppp distance cubes.
            The default is '0.1' [kpc]
        """
        self.path_to_bde = None
        self.path_to_file = None
        self.path_to_table = None
        self.path_to_input_table = None
        self.input_table = None
        self.save_input_table = False
        self.verbose = True
        self.gpy_setting = False
        self.intensity_threshold = 0.1
        self.distance_spacing = 0.1  # in [kpc]
        self.add_kinematic_distance = True
        self.colname_lon, self.colname_lat, self.colname_vel,\
            self.colname_kda = (None for i in range(4))
        self.colnr_lon, self.colnr_lat, self.colnr_vel,\
            self.colnr_kda = (None for i in range(4))
        self.prob_sa, self.prob_kd, self.prob_gl, self.prob_ps =\
            0.5, 1.0, 1.0, 0.25

        self.use_ncpus = None

    def set_probability_controls(self):
        s = '      '

        cwd = os.getcwd()
        os.chdir(self.path_to_bde)

        with open(os.path.join(
                self.path_to_bde, 'probability_controls.inp'), 'r') as fin:
            file_content = fin.readlines()
        with open(os.path.join(
                self.path_to_bde, 'probability_controls.inp'), 'w') as fout:
            for line in file_content:
                if not line.startswith('!'):
                    line = '{s}{a}{s}{b}{s}{c}{s}{d}'.format(
                        s=s, a=self.prob_sa, b=self.prob_kd, c=self.prob_gl,
                        d=self.prob_ps)
                fout.write(line)
        os.chdir(cwd)

    def make_fortran_out(self, source):
        """
        Create a fortran executable for the source.

        Replaces the default input file in the fortran script of the Bayesian
        distance estimator with the input file of the source, then creates a
        Fortran executable file.
        """
        with open("{}.f".format(self.path_to_source), "w") as fout:
            for line in self.bde_script:
                fout.write(line.replace('sources_info.inp',
                                        '{}_sources_info.inp'.format(source)))
        os.system('gfortran {}.f -o {}.out'.format(
                self.path_to_source, self.path_to_source))

    def extract_string(self, s, first, last, incl=False):
        """
        Search for a substring inside a string.

        Parameters
        ----------
        s : string that is searched for the substring
        first : first characters of the substring
        last : last characters of the substring
        incl : defines if the `first` and `last` characters are still part of
            the substring that will be returned. The default is `False`
            (`first` and `last` are not part of the returned substring)

        Returns
        -------
        substring of s
        """
        try:
            if incl is True:
                start = s.index(first)
                end = s.index(last) + len(last)
            else:
                start = s.index(first) + len(first)
                end = s.index(last, start)
            return s[start:end]
        except ValueError:
            return ""

    def extract_probability_info(self, line, lon, lat, pFar, kinDist):
        """
        Extract the distance results from the corresponding string in
        the output file of the Bayesian distance estimator tool.
        """
        deleteString = self.extract_string(
                line, 'Probability component', ':', incl=True)
        replaceString = self.extract_string(
                line, 'Probability component', ':')
        line = line.replace(deleteString, replaceString)
        line = line.replace('\n', '')
        comp, dist, err, prob, arm = line.split()
        c_u, c_v, c_w = self.get_cartesian_coords(lon, lat, float(dist))
        # if np.isnan(dist) is True:
        #     dist, err, prob = (0.0 for i in range(3))
        return [comp, dist, err, prob, arm,
                c_u, c_v, c_w, pFar, kinDist[0], kinDist[1]]

    def extract_results(self, result_file_content, kinDist=None):
        """
        Loop through the lines of the output file of the Bayesian distance
        estimator tool and search for the distance results.

        Parameters
        ----------
        result_file_content : List containing read-in lines of the output file
            ({source_name}.prt) of the Bayesian distance estimator tool
        """
        results = []
        flag = False
        for line in result_file_content:
            if flag:
                params = line.split()
                lon, lat, pFar =\
                    float(params[1]), float(params[2]), float(params[4])
                flag = False
            if 'Extra_info' in line:
                flag = True
            searchString = 'Probability component'
            if searchString in line:
                results.append(self.extract_probability_info(
                    line, lon, lat, pFar, kinDist))
        return results

    def extract_kinematic_distances(self, result_file_content):
        """"""
        kinDist = [np.NAN, np.NAN]

        flag = 'one'
        for line in result_file_content:
            searchString = 'Kinematic distance(s):'
            if searchString in line:
                if flag == 'one':
                    kinDist[0] = self.extract_kinematic_info(line)
                    flag = 'two'
                elif flag == 'two':
                    kinDist[1] = self.extract_kinematic_info(line)
        return kinDist

    def extract_kinematic_info(self, line):
        """
        Extract the distance results from the corresponding string in
        the output file of the Bayesian distance estimator tool.
        """
        line = line.replace('Kinematic distance(s):', '')
        line = line.replace('\n', '')
        return float(line)

    def get_results(self, source):
        """
        Extract the distance results from the output file ({source_name}.prt)
        of the Bayesian distance estimator tool.
        """
        for filename in [f for f in os.listdir(self.path_to_bde)
                         if f.startswith(source) and f.endswith(".prt")]:
            with open(os.path.join(self.path_to_bde, filename), 'r') as fin:
                result_file_content = fin.readlines()
        for filename in [f for f in os.listdir(self.path_to_bde) if f.startswith(source)]:
            os.remove(os.path.join(self.path_to_bde, filename))

        if self.add_kinematic_distance:
            kinDist = self.extract_kinematic_distances(result_file_content)
        else:
            kinDist = []
        results = self.extract_results(result_file_content, kinDist)
        return results

    def determine(self, row, idx):
        row = list(row)
        """
        Determine the distance of an lbv data point with the Bayesian distance
        estmator tool.
        """
        if self.gpy_setting:
            x_pos, y_pos, z_pos, intensity, lon, lat, vel = row
            source = "X{}Y{}Z{}".format(x_pos, y_pos, z_pos)
        else:
            # source, lon, lat, vel = row
            # source = "LON{}LAT{}VEL{}".format(
            #     row[self.colnr_lon], row[self.colnr_lat], row[self.colnr_vel])
            source = "SRC{}".format(str(idx).zfill(9))
            lon, lat, vel =\
                row[self.colnr_lon], row[self.colnr_lat], row[self.colnr_vel]

        if self.colnr_kda is not None:
            if row[self.colnr_kda] == 'F':
                p_far = 1.0
            elif row[self.colnr_kda] == 'N':
                p_far = 0.0
            else:
                p_far = 0.5
        else:
            p_far = self.check_KDA(lon, lat, vel)

        inputString = "{a}\t{b}\t{c}\t{d}\t{e}\t-\n".format(
            a=source, b=lon, c=lat, d=vel, e=p_far)
        self.path_to_source = os.path.join(self.path_to_bde, source)
        filepath = '{}_sources_info.inp'.format(self.path_to_source)
        with open(filepath, 'w') as fin:
            fin.write(inputString)

        self.make_fortran_out(source)
        cwd = os.getcwd()
        os.chdir(self.path_to_bde)
        os.system('{}.out'.format(source))
        os.chdir(cwd)

        rows = []
        if self.gpy_setting:
            row = [x_pos, y_pos, z_pos, intensity, lon, lat, vel]
        # else:
        #     row = [source, lon, lat, vel]
        results = self.get_results(source)
        for result in results:
            rows.append(row + result)

        return rows

    def check_KDA(self, lon, lat, vel):
        pFarVal = 0.5
        first = True
        found_entry = False
        for tableKdaInfo in [self.table_kda_info_3, self.table_kda_info_1, self.table_kda_info_2]:
            if not found_entry:
                for lonMin, lonMax, latMin, latMax, velMin, velMax, kda, pFar in zip(
                        tableKdaInfo['lonMin'], tableKdaInfo['lonMax'],
                        tableKdaInfo['latMin'], tableKdaInfo['latMax'],
                        tableKdaInfo['velMin'], tableKdaInfo['velMax'],
                        tableKdaInfo['KDA'], tableKdaInfo['pFar']):
                    if lonMin < lon < lonMax:
                        if latMin < lat < latMax:
                            if velMin < vel < velMax:
                                if first:
                                    kdaVal, pFarVal = kda, pFar
                                    first = False
                                else:
                                    if kdaVal != kda:
                                        pFarVal = 0.5

        return pFarVal

    def determine_column_indices(self):
        self.colnr_lon = self.input_table.colnames.index(self.colname_lon)
        self.colnr_lat = self.input_table.colnames.index(self.colname_lat)
        self.colnr_vel = self.input_table.colnames.index(self.colname_vel)
        if self.colname_kda is not None:
            self.colnr_kda = self.input_table.colnames.index(self.colname_kda)

    # def get_cartesian_coords(self, row):
    #     from astropy.coordinates import SkyCoord
    #     from astropy import units as u
    #
    #     c = SkyCoord(l=row[self.colname_lon]*u.degree,
    #                  b=row[self.colname_lat]*u.degree,
    #                  distance=row['dist']*u.kpc,
    #                  frame='galactic')
    #     c.representation = 'cartesian'
    #     c_u = round(c.u.value, 4)
    #     c_v = round(c.v.value, 4)
    #     c_w = round(c.w.value, 4)
    #
    #     return [c_u, c_v, c_w]

    def get_cartesian_coords(self, lon, lat, dist):
        from astropy.coordinates import SkyCoord
        from astropy import units as u

        c = SkyCoord(l=lon*u.degree,
                     b=lat*u.degree,
                     distance=dist*u.kpc,
                     frame='galactic')
        c.representation = 'cartesian'
        c_u = round(c.u.value, 4)
        c_v = round(c.v.value, 4)
        c_w = round(c.w.value, 4)

        return c_u, c_v, c_w

    def calculate_distances(self):
        self.check_settings()

        if self.verbose:
            string = str("prob_sa: {a}\nprob_kd: {b}\n"
                         "prob_gl: {c}\nprob_ps: {d}\n".format(
                             a=self.prob_sa, b=self.prob_kd, c=self.prob_gl,
                             d=self.prob_ps))
            print("setting probability controls to the following values:")
            print(string)

        self.set_probability_controls()

        if self.verbose:
            print('calculating Bayesian distance...')

        if self.gpy_setting:
            self.create_input_table()
        else:
            if self.input_table is None:
                self.input_table = Table.read(
                    self.path_to_input_table, format='ascii')
            self.determine_column_indices()

        import BD_wrapper.BD_multiprocessing as BD_multiprocessing
        BD_multiprocessing.init([self, self.input_table])
        results_list = BD_multiprocessing.func(use_ncpus=self.use_ncpus)
        print('SUCCESS\n')

        results_list = np.array([item for sublist in results_list
                                 for item in sublist])

        self.create_astropy_table(results_list)

    def initialize_data(self):
        self.dirname = os.path.dirname(self.path_to_file)
        self.file = os.path.basename(self.path_to_file)
        self.filename, self.fileExtension = os.path.splitext(self.file)

        self.dirname_table = os.path.dirname(self.path_to_table)
        self.table_file = os.path.basename(self.path_to_table)
        if not os.path.exists(self.dirname_table):
            os.makedirs(self.dirname_table)

        hdu = fits.open(self.path_to_file)[0]
        self.data = hdu.data
        self.header = hdu.header
        self.shape = (self.data.shape[0], self.data.shape[1],
                      self.data.shape[2])

    def check_settings(self):
        if self.path_to_bde is None:
            raise Exception("Need to specify 'path_to_bde'")
        path_to_file = os.path.join(
                self.path_to_bde, "Bayesian_distance_v1.0.f")
        with open(path_to_file, "r") as fin:
            bde_script = fin.readlines()
        self.bde_script = bde_script

        # if (self.path_to_file is None) and (self.path_to_input_table is None):
        #     errorMessage = str("specify 'path_to_file'")
        #     raise Exception(errorMessage)

        if self.path_to_table is None:
            errorMessage = str("specify 'path_to_table'")
            raise Exception(errorMessage)

        dirname = os.path.dirname(os.path.realpath(__file__))
        self.table_kda_info_1 = Table.read(
            os.path.join(dirname, 'KDA_info', 'KDA_info_EB+15.dat'),
            format='ascii')
        self.table_kda_info_2 = Table.read(
            os.path.join(dirname, 'KDA_info', 'KDA_info_RD+09.dat'),
            format='ascii')
        self.table_kda_info_3 = Table.read(
            os.path.join(dirname, 'KDA_info', 'KDA_info_Urquhart+17.dat'),
            format='ascii')

        self.dirname_table = os.path.dirname(self.path_to_table)
        if len(self.dirname_table) == 0:
            self.dirname_table = os.getcwd()
        self.table_file = os.path.basename(self.path_to_table)
        self.table_filename, self.table_file_extension =\
            os.path.splitext(self.table_file)
        if not os.path.exists(self.dirname_table):
            os.makedirs(self.dirname_table)

        text = 'Python wrapper for Bayesian distance estimator'
        border = len(text) * '='
        heading = '\n{a}\n{b}\n{a}\n'.format(a=border, b=text)
        if self.verbose:
            print(heading)

    def create_input_table(self):
        if self.verbose:
            print('creating input table...')

        self.initialize_data()

        velocityOffset = self.header['CRVAL3'] -\
            self.header['CDELT3']*(self.header['CRPIX3'] - 1)

        x_pos, y_pos, z_pos, intensity, longitude, latitude, velocity = (
                [] for i in range(7))

        for (x, y, z) in itertools.product(range(self.data.shape[2]),
                                           range(self.data.shape[1]),
                                           range(self.data.shape[0])):
            if float(self.data[z, y, x]) > self.intensity_threshold:
                x_pos.append(x)
                y_pos.append(y)
                z_pos.append(z)
                intensity.append(self.data[z, y, x])
                lon = (x - self.header['CRPIX1'])*self.header['CDELT1'] +\
                    self.header['CRVAL1']
                longitude.append(lon)
                lat = (y - self.header['CRPIX2'])*self.header['CDELT2'] +\
                    self.header['CRVAL2']
                latitude.append(lat)
                vel = (velocityOffset + np.array(z) *
                       self.header['CDELT3']) / 1000
                velocity.append(vel)

        names = ['x_pos', 'y_pos', 'z_pos', 'intensity', 'lon', 'lat', 'vel']
        self.input_table = Table([x_pos, y_pos, z_pos, intensity, longitude,
                                 latitude, velocity], names=names)

        if self.save_input_table:
            filename = '{}_input.dat'.format(self.filename)
            path_to_table = os.path.join(self.dirname_table, filename)
            self.input_table.write(path_to_table, format='ascii', overwrite=True)
            if self.verbose:
                print(">> saved input table '{}' in {}".format(
                        filename, self.dirname_table))

    def create_astropy_table(self, results):
        if self.verbose:
            print('creating Astropy table...')
        if self.gpy_setting:
            names = ('x_pos', 'y_pos', 'z_pos', 'intensity', 'lon', 'lat',
                     'vel', 'comp', 'dist', 'e_dist', 'prob', 'arm')
            dtype = ('i4', 'i4', 'i4', 'f4', 'f4', 'f4', 'f4',
                     'i4', 'f4', 'f4', 'f4', 'object')
        else:
            addedColnames = ['comp', 'dist', 'e_dist', 'prob', 'arm',
                             'c_u', 'c_v', 'c_w', 'pFar']
            if self.add_kinematic_distance:
                addedColnames += ['kDist_1', 'kDist_2']
            names = self.input_table.colnames + addedColnames

            dtypeinput_table = []
            for name, dtype in self.input_table.dtype.descr:
                dtypeinput_table.append(dtype)
            added_dtype = ['i4', 'f4', 'f4', 'f4', 'object',
                           'f4', 'f4', 'f4', 'f4']
            if self.add_kinematic_distance:
                added_dtype += ['f4', 'f4']
            dtype = dtypeinput_table + added_dtype

        self.table_results = Table(data=results, names=names, dtype=dtype)

        for key in ['dist', 'e_dist', 'prob', 'c_u', 'c_v', 'c_w']:
            if key in self.table_results.colnames:
                self.table_results[key].format = "{0:.4f}"
        for key in ['pFar', 'kDist_1', 'kDist_2']:
            if key in self.table_results.colnames:
                self.table_results[key].format = "{0:.2f}"

        if self.verbose:
            print(">> saved table '{}' in {}\n".format(
                    self.table_file, self.dirname_table))

        self.table_results.write(self.path_to_table, format='ascii',
                                overwrite=True)

    def get_table_distance_max_probability(self):
        from tqdm import tqdm
        if self.verbose:
            print('creating Astropy table containing only distance result '
                  'with the highest probability...')

        # remove_rows = []
        #
        # for idx, component in tqdm(enumerate(self.tableResults['comp'])):
        #     if idx == 0:
        #         comps_indices = [idx]
        #     else:
        #         if component == 1:
        #             first_idx = True
        #             for comps_idx in comps_indices:
        #                 if first_idx:
        #                     prob = self.tableResults['prob'][comps_idx]
        #                     prob_idx = comps_idx
        #                     first_idx = False
        #                 else:
        #                     if self.tableResults['prob'][comps_idx] < prob:
        #                         remove_rows.append(comps_idx)
        #                     else:
        #                         remove_rows.append(prob_idx)
        #                         prob_idx = comps_idx
        #                         prob = self.tableResults['prob'][comps_idx]
        #             comps_indices = [idx]
        #         if component != 1:
        #             comps_indices.append(idx)

        remove_rows = np.array([])

        for idx, component in tqdm(enumerate(self.table_results['comp'])):
            if idx == 0:
                comps_indices = np.array([idx])
            else:
                if (component == 1):
                    if comps_indices.size > 1:
                        sort_indices_highest_probability = np.argsort(
                            self.table_results['prob'][comps_indices])[::-1]
                        remove = sort_indices_highest_probability[1:]
                        remove_rows = np.append(remove_rows, comps_indices[remove])
                    comps_indices = np.array([idx])
                else:
                    comps_indices = np.append(comps_indices, idx)

        #  take care of the last distance results in the list
        sort_indices_highest_probability = np.argsort(
            self.table_results['prob'][comps_indices])[::-1]
        remove = sort_indices_highest_probability[1:]
        remove_rows = np.append(remove_rows, comps_indices[remove])

        remove_rows = remove_rows.astype(int)
        self.table_results.remove_rows(remove_rows)

        self.table_file = '{}{}{}'.format(self.table_filename, '_p_max',
                                          self.table_file_extension)
        self.path_to_table = os.path.join(self.dirname_table, self.table_file)

        if self.verbose:
            print(">> saved table '{}' in {}".format(
                    self.table_file, self.dirname_table))

        self.table_results.write(self.path_to_table, format='ascii',
                                 overwrite=True)

    def find_index_max_probability(self, indices, arm=False):
        idx = [i for i in indices]
        prob = [self.table['prob'][i] for i in indices]
        max_idx = prob.index(max(prob))

        if arm:
            arms = [self.table['arm'][i] for i in indices]
            return idx[max_idx], arms[max_idx]
        else:
            return idx[max_idx]

    def make_ppp_intensity_cube(self):
        if self.verbose:
            print('create PPP weighted intensity cube...')

        self.check_settings()
        self.initialize_data()

        self.table = Table.read(self.path_to_table, format='ascii.fixed_width')
        maxDist = int(max(self.table['dist'])) + 1
        zrange = int(maxDist/self.distance_spacing)
        self.shape = (zrange, self.data.shape[1], self.data.shape[2])
        array = np.zeros(self.shape, dtype='float32')
        self.header['NAXIS3'] = zrange
        self.header['CRPIX3'] = 1.
        self.header['CRVAL3'] = self.distance_spacing
        self.header['CDELT3'] = self.distance_spacing
        self.header['CTYPE3'] = 'DISTANCE'
        index_list = []

        for idx, (component, probability) in enumerate(
                zip(self.table['comp'], self.table['dist'])):
            if idx == 0:
                comps_indices = [idx]
            else:
                if component == 1:
                    index = self.find_index_max_probability(comps_indices)
                    index_list.append(index)

                    x = self.table['x_pos'][index]
                    y = self.table['y_pos'][index]
                    dist = round(self.table['dist'][index], 1)
                    z = round(dist / self.distance_spacing)
                    intensity = self.table['intensity'][index]

                    array[z, y, x] += intensity

                    comps_indices = [idx]
                if component != 1:
                    comps_indices.append(idx)

        filename = '{}_distance_ppp.fits'.format(self.filename)
        pathname = os.path.join(self.dirname_table, 'FITS')
        if not os.path.exists(pathname):
            os.makedirs(pathname)
        path_to_file = os.path.join(pathname, filename)
        fits.writeto(path_to_file, array, self.header, overwrite=True)
        if self.verbose:
            print(">> saved '{}' to {}".format(filename, pathname))

        def make_ppv_distance_cube(self):
            if self.verbose:
                print('create PPV distance cube...')

            self.check_settings()
            self.initialize_data()

            self.table = Table.read(self.path_to_table,
                                    format='ascii.fixed_width')
            array = np.zeros(self.shape, dtype='float32')
            index_list = []

            for idx, (component, probability) in enumerate(
                    zip(self.table['comp'], self.table['dist'])):
                if idx == 0:
                    comps_indices = [idx]
                else:
                    if component == 1:
                        index = self.find_index_max_probability(comps_indices)
                        index_list.append(index)

                        x = self.table['x_pos'][index]
                        y = self.table['y_pos'][index]
                        z = self.table['z_pos'][index]
                        dist = self.table['dist'][index]

                        array[z, y, x] = dist

                        comps_indices = [idx]
                    if component != 1:
                        comps_indices.append(idx)

            filename = '{}_distance.fits'.format(self.filename)
            pathname = os.path.join(self.dirname_table, 'FITS')
            if not os.path.exists(pathname):
                os.makedirs(pathname)
            path_to_file = os.path.join(pathname, filename)
            fits.writeto(path_to_file, array, self.header, overwrite=True)
            if self.verbose:
                print(">> saved '{}' to {}".format(filename, pathname))

    def timer(self, mode='start', start_time=None):
        """"""
        import time

        if mode == 'start':
            return time.time()
        elif mode == 'stop':
            print('\njob finished on {}'.format(time.ctime()))
            print('required run time: {:.4f} s\n'.format(
                time.time() - start_time))
