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
        self.path_to_input_table = None
        self.path_to_output_table = None
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
        self.plot_probability = False

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
                             kin_dist=None, kda_ref=None):
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

                if kda_ref is not None:
                    result += [kda_ref]

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
                             kin_dist=None, kda_ref=None):
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

                if kda_ref is not None:
                    result += [kda_ref]

                if kin_dist is not None:
                    result += kin_dist

                results.append(result)
        return results

    def get_results(self, source, kda_ref):
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

        results = self._p[self.version]['fct_extract'](
            input_file_content, result_file_content, kin_dist=kinDist, kda_ref=kda_ref)

        if self.plot_probability:
            self.plot_probability_density(source, results)

        #  remove for TESTING:
        if not self.save_temporary_files:
            for filename in [f for f in os.listdir(self.path_to_bde) if f.startswith(source)]:
                os.remove(os.path.join(self.path_to_bde, filename))

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
        os.system('./{}.out'.format(source))
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
        kda_ref = None

        if self.colnr_kda is not None:
            if row[self.colnr_kda] == 'F':
                p_far = 1.0
            elif row[self.colnr_kda] == 'N':
                p_far = 0.0
            elif isinstance(row[self.colnr_kda], float):
                p_far = row[self.colnr_kda]
        elif self.check_for_kda_solutions:
            p_far, kda_ref = self.check_KDA(lon, lat, vel)

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
        results = self.get_results(source, kda_ref)
        for result in results:
            rows.append(row + result)

        return rows

    def get_values_from_init_file(self, init_file):
        """Read in values from init file.

        """
        import ast
        import configparser
        config = configparser.ConfigParser()
        config.read(init_file)

        for key, value in config['DEFAULT'].items():
            try:
                setattr(self, '_' + key, ast.literal_eval(value))
            except ValueError:
                raise Exception('Could not parse parameter {} from config file'.format(key))

    def point_in_ellipse(self, table, lon, lat):
        cos_pa = table['cos_pa'].data.data
        sin_pa = table['sin_pa'].data.data
        glon = table['GLON'].data.data
        glat = table['GLAT'].data.data
        aa = table['aa'].data.data
        bb = table['bb'].data.data

        a = (cos_pa * (lon - glon) + sin_pa * (lat - glat))**2
        b = (sin_pa * (lon - glon) - cos_pa * (lat - glat))**2
        ellipse = (a / aa) + (b / bb)
        ellipse[ellipse < 1] = 1
        weight = 1 / ellipse

        return weight >= self._threshold_spatial, weight

    def get_weight_velocity(self, table, vel):
        """Calculate the weight for the velocity association.

        Gaussian function: amp * np.exp(-4. * np.log(2) * (x-mean)**2 / fwhm**2)
        mean = 0
        Conversion factor std to FWHM = 2.354820045

        Renormalization factor for amplitude, so that Gaussian function is 1 at the standard deviation; scale = 1 / (np.exp(-4. * np.log(2) / (2.354820045)**2)) = 1.6487212707217973

        """
        normalized_amp = 1.6487212707217973
        fwhm_factor = 2.354820045

        vlsr = table['VLSR'].data.data
        dvlsr = table['d_VLSR'].data.data

        x = np.abs(vlsr - vel) / dvlsr

        weight = normalized_amp * np.exp(
            -4. * np.log(2) * x**2 / fwhm_factor**2)
        weight[weight > 1] = 1

        return weight >= self._threshold_spectral, weight

    def check_KDA(self, lon, lat, vel):
        p_far = 0.5
        ref = '--'
        dirname = os.path.dirname(os.path.realpath(__file__))
        for table, tablename in zip(self._kda_tables, self.kda_info_tables):
            path_to_ini = os.path.join(dirname, 'KDA_info', tablename + '.ini')
            self.get_values_from_init_file(path_to_ini)
            mask_pp, weight_pp = self.point_in_ellipse(table, lon, lat)
            mask_vlsr, weight_vlsr = self.get_weight_velocity(table, vel)

            mask_total = np.logical_and(mask_pp, mask_vlsr)
            weight_total = weight_pp * weight_vlsr
            weight_total = weight_total[mask_total]
            p_far_values = table['p_far'].data
            p_far_values = p_far_values[mask_total]

            n_values = np.count_nonzero(mask_total)
            if n_values == 0:
                continue
            elif n_values == 1:
                p_far = 0.5 + self._weight_pfar * p_far_values * weight_total
                print('single:', p_far)
                ref = self._reference
                break
            else:
                p_far = 0.5 + self._weight_pfar * (np.average(
                    p_far_values * weight_total, weights=weight_total))
                print('multiple:', p_far)
                ref = self._reference
                break

        return round(float(p_far), 2), ref

    def determine_column_indices(self):
        self.colnr_lon = self.input_table.colnames.index(self.colname_lon)
        self.colnr_lat = self.input_table.colnames.index(self.colname_lat)
        self.colnr_vel = self.input_table.colnames.index(self.colname_vel)
        if self.colname_e_vel is not None:
            self.colnr_e_vel = self.input_table.colnames.index(self.colname_e_vel)
        if self.colname_kda is not None:
            self.colnr_kda = self.input_table.colnames.index(self.colname_kda)

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

        if self.path_to_output_table is not None:
            self.path_to_table = self.path_to_output_table

        if self.path_to_table is None:
            errorMessage = str("specify 'path_to_table'")
            raise Exception(errorMessage)

        dirname = os.path.dirname(os.path.realpath(__file__))
        self._kda_tables, self._kda_tables_ref = [], []
        for tablename in self.kda_info_tables:
            self._kda_tables.append(
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
            added_colnames = ['comp', 'dist', 'e_dist', 'prob', 'arm',
                              'c_u', 'c_v', 'c_w', 'p_far']

            dtypeinput_table = []
            for name, dtype in self.input_table.dtype.descr:
                dtypeinput_table.append(dtype)
            added_dtype = ['i4', 'f4', 'f4', 'f4', 'object',
                           'f4', 'f4', 'f4', 'f4']

            if self.check_for_kda_solutions and (self.colname_kda is None):
                added_colnames += ['KDA_ref']
                added_dtype += ['object']
            if self.add_kinematic_distance:
                added_colnames += ['kDist_1', 'kDist_2']
                added_dtype += ['f4', 'f4']
            names = self.input_table.colnames + added_colnames
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

    def plot_probability_density(self, source, results):
        import matplotlib.pyplot as plt

        def get_maximum_distance(distance, probability, max_dist=None):
            try:
                prob_threshold = 0.05
                distance_cutoff = distance[probability > prob_threshold][-1]
                return max(max_dist, int(distance_cutoff + 1))
            except IndexError:
                return max_dist

        distKD, probKD = np.loadtxt(
            os.path.join(self.path_to_bde, '{}_kinematic_distance_pdf.dat'.format(source)),
            usecols=(0, 1), skiprows=2, unpack=True)

        if self.version == '1.0':
            distGL, probGL = np.loadtxt(
                os.path.join(self.path_to_bde, '{}_latitude_pdf.dat'.format(source)),
                usecols=(0, 1), skiprows=2, unpack=True)

            distSA, probSA = np.loadtxt(
                os.path.join(self.path_to_bde, '{}_spiral_arm_pdf.dat'.format(source)),
                usecols=(0, 1), skiprows=2, unpack=True)
        elif self.version == '2.4':
            distSA, probSA = np.loadtxt(
                os.path.join(self.path_to_bde, '{}_arm_latitude_pdf.dat'.format(source)),
                usecols=(0, 1), skiprows=2, unpack=True)

            distGL, probGL = distSA, probSA

        arm_ranges_lower, arm_ranges_upper = np.loadtxt(
            os.path.join(self.path_to_bde, '{}_arm_ranges.dat'.format(source)),
            usecols=(0, 1), skiprows=2, unpack=True)

        spiral_arms = np.genfromtxt(
            os.path.join(self.path_to_bde, '{}_arm_ranges.dat'.format(source)),
            skip_header=2, usecols=2, dtype='str')

        skip_spiral_arm_ranges = False
        if spiral_arms.size == 2:
            arm_ranges_lower = [arm_ranges_lower[1]]
            arm_ranges_upper = [arm_ranges_upper[1]]
            spiral_arms = [spiral_arms[1]]
        elif spiral_arms.size > 2:
            arm_ranges_lower = arm_ranges_lower[1:]
            arm_ranges_upper = arm_ranges_upper[1:]
            spiral_arms = spiral_arms[1:]
        else:
            skip_spiral_arm_ranges = True

        distPS, probPS = np.loadtxt(
            os.path.join(self.path_to_bde, '{}_parallaxes_pdf.dat'.format(source)),
            usecols=(0, 1), skiprows=2, unpack=True)

        distFD, probFD = np.loadtxt(
            os.path.join(self.path_to_bde, '{}_final_distance_pdf.dat'.format(source)),
            usecols=(0, 1), skiprows=2, unpack=True)

        max_dist = 0
        for dist, prob in zip(
                [distKD, distGL, distSA, distPS], [probKD, probGL, probSA, probPS]):
            max_dist = get_maximum_distance(dist, prob, max_dist=max_dist)

        fig = plt.figure(figsize=(10, 7.5))
        ax = fig.add_subplot(1, 1, 1)

        ax.set_xlabel('Distance [kpc]', size=20)
        ax.set_ylabel('Probability density [kpc$^{-1}$]', size=20)

        ax.tick_params(axis='both', labelsize=16, pad=8)
        ax.tick_params(axis='both', which='major', direction='out',
                       width=1.25, length=10, pad=8)
        ax.tick_params(axis='both', which='minor', direction='out',
                       width=1.25, length=5)

        ax.plot(distKD, probKD, label='KD', lw=2.5, ls='solid', c='dodgerblue', alpha=0.75)
        ax.plot(distGL, probGL, label='GL', lw=2.5, ls='dotted', c='orange', alpha=1.0)
        ax.plot(distSA, probSA, label='SA', lw=2.5, ls='--', c='indianred', alpha=1.0)
        ax.plot(distPS, probPS, label='PS', lw=2.5, ls='-.', c='forestgreen', alpha=1.0)
        ax.plot(distFD, probFD, label='combined', lw=2, ls='solid', c='black', alpha=1.0)

        if not skip_spiral_arm_ranges:
            for lower, upper, text in zip(
                    arm_ranges_lower, arm_ranges_upper, spiral_arms):
                if lower > max_dist:
                    continue
                ax.axvspan(lower, upper, alpha=0.15, color='indianred')
                ax.text((upper + lower)/2, ax.get_ylim()[1] * 0.99, text, size=16, color='indianred',
                        horizontalalignment='center', verticalalignment='top')

        box = ax.get_position()
        ax.set_position([box.x0, box.y0,
                         box.width, box.height * 0.9])
        # Put a legend below current axis
        leg1 = ax.legend(
            loc='upper center', bbox_to_anchor=(0.5, 1.08),
            fancybox=False, shadow=False, ncol=5,
            fontsize=14, numpoints=1, frameon=0)

        markers, texts = [], []

        for result in results:
            dist = float(result[1])
            if dist <= 0:
                continue
            e_dist = float(result[2])
            prob = float(result[3])
            text = '{a:.2f}$\\pm${b:.2f}kpc ({c:.0%})'.format(
                a=dist, b=e_dist, c=prob)
            marker = ax.scatter(dist, 0)
            markers.append(marker)
            texts.append(text)
            ax.errorbar(dist, 0, xerr=e_dist)

        leg2 = ax.legend(
            markers, texts,
            loc='upper center', bbox_to_anchor=(0.5, 1.135),
            fancybox=False, shadow=False, ncol=len(results),
            fontsize=14, numpoints=1, frameon=0)

        ax.set_xlim([0, max_dist])

        ax.add_artist(leg1)

        path_to_file = os.path.join(self.dirname_table, source + '.pdf')
        plt.savefig(path_to_file, bbox_inches='tight')
        plt.close()
