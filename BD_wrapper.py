import itertools
import os
import pickle

import numpy as np

from astropy.io import fits
from astropy.table import Table

from functools import reduce


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
        self.version = '2.4'
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
        self.check_for_kda_solutions = True
        self.colname_lon, self.colname_lat, self.colname_vel,\
            self.colname_e_vel, self.colname_kda = (None for i in range(5))
        self.colnr_lon, self.colnr_lat, self.colnr_vel,\
            self.colnr_e_vel, self.colnr_kda = (None for i in range(5))
        self.prob_sa, self.prob_kd, self.prob_gl, self.prob_ps, self.prob_pm =\
            (None for _ in range(5))
        self.table_format = 'ascii'
        self.save_temporary_files = False
        self.default_e_vel = 5.0
        self.kda_info_tables = ['Urquhart+18', 'Ellsworth-Bowers+15',
                                'Roman-Duval+09']

        self.use_ncpus = None

        self._p = {
            '1.0': {
                'bdc_fortran': 'Bayesian_distance_v1.0.f',
                'summary_suffix': '.prt',
                'fct_extract': self.extract_results_v1p0},
            '2.4': {
                'bdc_fortran': 'Bayesian_distance_2019_fromlist_v2.4.f',
                'summary_suffix': 'summary.prt',
                'fct_extract': self.extract_results_v2p4}
        }

    def say(self, message, end=None):
        """Diagnostic messages."""
        if self.verbose:
            print(message, end=end)

    def set_probability_controls(self):
        s = '      '

        default_vals = {
            '1.0': {'SA': 0.5, 'KD': 1.0, 'GL': 1.0, 'PS': 0.25, 'PM': None},
            '2.4': {'SA': 0.85, 'KD': 0.85, 'GL': 0.85, 'PS': 0.15, 'PM': 0.85}
            }

        if self.prob_sa is None:
            self.prob_sa = default_vals[self.version]['SA']
        if self.prob_kd is None:
            self.prob_kd = default_vals[self.version]['KD']
        if self.prob_gl is None:
            self.prob_gl = default_vals[self.version]['GL']
        if self.prob_ps is None:
            self.prob_ps = default_vals[self.version]['PS']
        if self.prob_pm is None:
            self.prob_pm = default_vals[self.version]['PM']

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
                    if self.prob_pm is not None:
                        line += '{s}{a}'.format(s=s, a=self.prob_pm)
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

    def extract_probability_info(self, line, lon, lat, p_far):
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
        return [comp, dist, err, prob, arm, c_u, c_v, c_w, p_far]

    def extract_results_v1p0(self, input_file_content, result_file_content,
                             kin_dist=None):
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
                lon, lat, p_far =\
                    float(params[1]), float(params[2]), float(params[4])
                flag = False
            if 'Extra_info' in line:
                flag = True
            searchString = 'Probability component'
            if searchString in line:
                result = self.extract_probability_info(line, lon, lat, p_far)

                if kin_dist is not None:
                    result += kin_dist

                results.append(result)
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

    def extract_results_v2p4(self, input_file_content, result_file_content,
                             kin_dist=None):
        for line in input_file_content:
            if line.startswith('!'):
                continue
            params = line.split()
            p_far = params[5]

        for line in result_file_content:
            if line.startswith('!'):
                continue
            params = line.split()

            n_params = len(params)

            lon, lat, vlsr, e_vlsr = params[:4]

            results = []

            for i in range(1, int(n_params / 4)):
                comp = int(n_params / 4) - 1
                dist, e_dist, prob, arm = params[i*4:(i + 1)*4]
                c_u, c_v, c_w = self.get_cartesian_coords(
                    float(lon), float(lat), float(dist))

                result = [comp, dist, e_dist, prob, arm, c_u, c_v, c_w, p_far]

                if kin_dist is not None:
                    result += kin_dist

                results.append(result)
        return results

    def get_results(self, source):
        """
        Extract the distance results from the output file ({source_name}.prt)
        of the Bayesian distance estimator tool.
        """
        suffix = self._p[self.version]['summary_suffix']
        for filename in [f for f in os.listdir(self.path_to_bde)
                         if f.startswith(source) and f.endswith(suffix)]:
            with open(os.path.join(self.path_to_bde, filename), 'r') as fin:
                result_file_content = fin.readlines()

        for filename in [f for f in os.listdir(self.path_to_bde)
                         if f.startswith(source) and f.endswith("info.inp")]:
            with open(os.path.join(self.path_to_bde, filename), 'r') as fin:
                input_file_content = fin.readlines()

        if self.add_kinematic_distance:
            if self.version == '1.0':
                kd_content = result_file_content.copy()
            elif self.version == '2.4':
                with open(os.path.join(self.path_to_bde, source + '.prt'), 'r') as fin:
                    kd_content = fin.readlines()
            kinDist = self.extract_kinematic_distances(kd_content)
        else:
            kinDist = None

        #  TESTING:
        for filename in [f for f in os.listdir(self.path_to_bde) if f.startswith(source)]:
            os.remove(os.path.join(self.path_to_bde, filename))

        results = self._p[self.version]['fct_extract'](
            input_file_content, result_file_content, kin_dist=kinDist)

        # if self.version == '1.0':
        #     results = self.extract_results_v1p0(result_file_content, kinDist)
        # elif self.version == '2.4':
        #     results = self.extract_results_v2p4(
        #         input_file_content, result_file_content)
        return results

    def run_bdc_script(self, source, input_string):
        self.path_to_source = os.path.join(self.path_to_bde, source)
        filepath = '{}_sources_info.inp'.format(self.path_to_source)
        with open(filepath, 'w') as fin:
            fin.write(input_string)
        self.make_fortran_out(source)
        cwd = os.getcwd()
        os.chdir(self.path_to_bde)
        os.system('{}.out'.format(source))
        os.chdir(cwd)

    def bdc_calculation_ok(self, source):
        suffix = self._p[self.version]['summary_suffix']
        for filename in [f for f in os.listdir(self.path_to_bde)
                         if f.startswith(source) and f.endswith(suffix)]:
            with open(os.path.join(self.path_to_bde, filename), 'r') as fin:
                result_file_content = fin.readlines()
        for line in result_file_content:
            if not line.startswith('!'):
                return True

        return False

    def determine(self, row, idx):
        row = list(row)
        """
        Determine the distance of an lbv data point with the Bayesian distance
        estmator tool.
        """
        e_vel = None

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

            if self.colnr_e_vel is not None:
                e_vel = row[self.colnr_e_vel]
                if abs(float(e_vel)) > 10:#abs(float(vel)):
                    e_vel = None

        p_far = 0.5

        if self.colnr_kda is not None:
            if row[self.colnr_kda] == 'F':
                p_far = 1.0
            elif row[self.colnr_kda] == 'N':
                p_far = 0.0
        elif self.check_for_kda_solutions:
            p_far = self.check_KDA(lon, lat, vel)

        if self.version == '1.0':
            plusminus = ''
        elif self.version == '2.4':
            # TODO: implement minimum error for velocity
            if e_vel is not None:
                plusminus = '{}\t'.format(e_vel)
            else:
                plusminus = '{}\t'.format(self.default_e_vel)

        input_string = "{a}\t{b}\t{c}\t{d}\t{e}{f}\t-\n".format(
            a=source, b=lon, c=lat, d=vel, e=plusminus, f=p_far)

        self.run_bdc_script(source, input_string)

        if (self.version == '2.4') and (p_far != 0.5):
            if not self.bdc_calculation_ok(source):
                for filename in [f for f in os.listdir(self.path_to_bde) if f.startswith(source)]:
                    os.remove(os.path.join(self.path_to_bde, filename))

                p_far = 0.5
                input_string = "{a}\t{b}\t{c}\t{d}\t{e}{f}\t-\n".format(
                    a=source, b=lon, c=lat, d=vel, e=plusminus, f=p_far)
                self.run_bdc_script(source, input_string)

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
        p_far = 0.5
        for table in self.kda_tables:
            lon_min, lon_max = table['lonMin'].data, table['lonMax'].data
            lat_min, lat_max = table['latMin'].data, table['latMax'].data
            vel_min, vel_max = table['velMin'].data, table['velMax'].data
            # kda, pFar = table['KDA'].data, table['pFar'].data

            i_lon = np.intersect1d(np.where(lon_min < lon)[0], np.where(lon_max > lon)[0])
            i_lat = np.intersect1d(np.where(lat_min < lat)[0], np.where(lat_max > lat)[0])
            i_vel = np.intersect1d(np.where(vel_min < vel)[0], np.where(vel_max > vel)[0])
            indices = reduce(np.intersect1d, (i_lon, i_lat, i_vel))

            n_values = len(indices)
            if n_values == 0:
                continue
            elif n_values == 1:
                p_far = table['pFar'][indices].data
                break
            else:
                p_far = sum(table['pFar'][indices].data) / n_values
                break

        return round(float(p_far), 2)

    def determine_column_indices(self):
        self.colnr_lon = self.input_table.colnames.index(self.colname_lon)
        self.colnr_lat = self.input_table.colnames.index(self.colname_lat)
        self.colnr_vel = self.input_table.colnames.index(self.colname_vel)
        if self.colname_e_vel is not None:
            self.colnr_e_vel = self.input_table.colnames.index(self.colname_e_vel)
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

        self.set_probability_controls()

        string = str("prob_sa: {a}\nprob_kd: {b}\n"
                     "prob_gl: {c}\nprob_ps: {d}\n".format(
                         a=self.prob_sa, b=self.prob_kd, c=self.prob_gl,
                         d=self.prob_ps))
        if self.version == '2.4':
            string += 'prob_pm: {}\n'.format(self.prob_pm)
        self.say("setting probability controls to the following values:")
        self.say(string)

        self.say('calculating Bayesian distance...')

        if self.gpy_setting:
            self.create_input_table()
        else:
            if self.input_table is None:
                self.input_table = Table.read(
                    self.path_to_input_table, format='ascii')
                #  TESTING:
                # self.input_table = self.input_table[62000:62001]
            self.determine_column_indices()

        import BD_wrapper.BD_multiprocessing as BD_multiprocessing
        BD_multiprocessing.init([self, self.input_table])
        results_list = BD_multiprocessing.func(use_ncpus=self.use_ncpus)
        print('SUCCESS\n')

        for i, item in enumerate(results_list):
            if not isinstance(item, list):
                self.say("Error for distance with index {}: {}".format(i, item))
                del results_list[i]
                continue

        results_list = np.array([item for sublist in results_list
                                 for item in sublist])

        if self.save_temporary_files:
            filepath = os.path.join(
                os.path.dirname(self.path_to_input_table),
                '_bdc_results_list.pickle')
            with open(filepath, 'wb') as p_file:
                pickle.dump(results_list, p_file)

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
        if (self.path_to_bde is None) and (self.version is None):
            raise Exception("Need to specify 'path_to_bde' or 'version'")

        path_script = os.path.dirname(os.path.realpath(__file__))

        if self.version is None:
            path_to_file = os.path.join(
                    self.path_to_bde, "Bayesian_distance_v1.0.f")
            if not os.path.exists(path_to_file):
                path_to_file = os.path.join(
                        self.path_to_bde, "Bayesian_distance_2019_fromlist_v2.4.f")
            self.version = self.extract_string(
                os.path.basename(path_to_file), '_v', '.f')
        else:
            self.path_to_bde = os.path.join(
                path_script, 'BDC', 'v' + self.version)
            path_to_file = os.path.join(
                self.path_to_bde, self._p[self.version]['bdc_fortran'])

        with open(path_to_file, "r") as fin:
            bde_script = fin.readlines()
        self.bde_script = bde_script

        if self.path_to_table is None:
            errorMessage = str("specify 'path_to_table'")
            raise Exception(errorMessage)

        dirname = os.path.dirname(os.path.realpath(__file__))
        self.kda_tables = []
        for tablename in self.kda_info_tables:
            self.kda_tables.append(
                Table.read(
                    os.path.join(dirname, 'KDA_info', tablename + '.dat'),
                    format='ascii')
            )

        self.dirname_table = os.path.dirname(self.path_to_table)
        if len(self.dirname_table) == 0:
            self.dirname_table = os.getcwd()
        self.table_file = os.path.basename(self.path_to_table)
        self.table_filename, self.table_file_extension =\
            os.path.splitext(self.table_file)
        if not os.path.exists(self.dirname_table):
            os.makedirs(self.dirname_table)

        text = 'Python wrapper for Bayesian distance calculator v{}'.format(
            self.version)
        border = len(text) * '='
        heading = '\n{a}\n{b}\n{a}\n'.format(a=border, b=text)
        self.say(heading)

    def create_input_table(self):
        self.say('creating input table...')

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
            self.say(">> saved input table '{}' in {}".format(
                     filename, self.dirname_table))

    def create_astropy_table(self, results):
        self.say('creating Astropy table...')
        if self.gpy_setting:
            names = ('x_pos', 'y_pos', 'z_pos', 'intensity', 'lon', 'lat',
                     'vel', 'comp', 'dist', 'e_dist', 'prob', 'arm')
            dtype = ('i4', 'i4', 'i4', 'f4', 'f4', 'f4', 'f4',
                     'i4', 'f4', 'f4', 'f4', 'object')
        else:
            addedColnames = ['comp', 'dist', 'e_dist', 'prob', 'arm',
                             'c_u', 'c_v', 'c_w', 'p_far']
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
        for key in ['p_far', 'kDist_1', 'kDist_2']:
            if key in self.table_results.colnames:
                self.table_results[key].format = "{0:.2f}"

        self.say(">> saved table '{}' in {}\n".format(
                 self.table_file, self.dirname_table))

        self.table_results.write(self.path_to_table, format=self.table_format,
                                 overwrite=True)

    def get_table_distance_max_probability(self, save=True):
        from tqdm import tqdm
        self.say('creating Astropy table containing only distance results '
                 'with the highest probability...')

        remove_rows = np.array([])

        if self.version == '1.0':
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
        elif self.version == '2.4':
            comps_indices = np.array([], dtype='int')

            for idx, component in tqdm(enumerate(self.table_results['comp'])):
                comps_indices = np.append(comps_indices, idx)

                if len(comps_indices) == component:
                    #  TODO: in case of 50/50 split of components the first one gets discarded by default!
                    remove = np.argmin(self.table_results['prob'][comps_indices])
                    remove_rows = np.append(remove_rows, comps_indices[remove])
                    comps_indices = np.array([], dtype='int')

        remove_rows = remove_rows.astype(int)
        self.table_results.remove_rows(remove_rows)

        if save:
            self.table_file = '{}{}{}'.format(self.table_filename, '_p_max',
                                              self.table_file_extension)
            self.path_to_table = os.path.join(
                self.dirname_table, self.table_file)

            self.say(">> saved table '{}' in {}".format(
                     self.table_file, self.dirname_table))

            self.table_results.write(self.path_to_table,
                                     format=self.table_format,
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
        self.say('create PPP weighted intensity cube...')

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
        self.say(">> saved '{}' to {}".format(filename, pathname))

        def make_ppv_distance_cube(self):
            self.say('create PPV distance cube...')

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
            self.say(">> saved '{}' to {}".format(filename, pathname))

    def timer(self, mode='start', start_time=None):
        """"""
        import time

        if mode == 'start':
            return time.time()
        elif mode == 'stop':
            print('\njob finished on {}'.format(time.ctime()))
            print('required run time: {:.4f} s\n'.format(
                time.time() - start_time))
