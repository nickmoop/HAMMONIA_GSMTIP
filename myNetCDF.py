from datetime import datetime, timedelta
import os

from netCDF4 import Dataset


PRESSURE = 101325 * 2.7 ** (80000 / -7000)
GAS_CONSTANT = 8.3144598


def process_file_ik19(filename):
    base_dir = "/home/nick/python/netCDF"
    file_to_read = open(filename, 'r')
    count = 0
    tmp = 0
    all_lines = file_to_read.readlines()
    i = 0
    current_date = 0
    while i < len(all_lines):
        tmp += 1
        line = all_lines[i]
        if 'seans' in line or 'SEANS' in line:
            date_str = all_lines[i + 1].split(',')[0].replace('\n', '').replace('"', '')
            if len(date_str.split('.')[2]) == 4:
                current_date = datetime.datetime.strptime(date_str, '%d.%m.%Y').date()
            else:
                current_date = datetime.datetime.strptime(date_str, '%d.%m.%y').date()
            print(current_date)
            count += 1
            i += 1
        else:
            data_str = line.replace('\n', '').split(',')
            altitude = 0
            frequency = 0
            density = 0
            #print(data_str)
            if len(data_str) == 10 and data_str[9] != 'NaN':
                frequency = float(data_str[9])
                density = float(data_str[9])
            if len(data_str) == 11 and data_str[10] != 'NaN':
                altitude = float(data_str[10])
            if altitude != 0 or frequency != 0:
                current_time = datetime.datetime.strptime('{:06d}'.format(int(data_str[2])), '%H%M%S').time()
                filename = make_filename(current_date.year, current_date.timetuple().tm_yday, current_time.hour, current_time.minute, data_str[0], data_str[1])
                latitude = float(data_str[5])
                longitude = float(data_str[6])
                #print(altitude, frequency)
                path_to_write = '{}/IK19netCDF/{}/{}.{:03d}/'.format(base_dir, current_date.year, current_date.year, current_date.timetuple().tm_yday)
                if not os.path.exists(path_to_write):
                    os.makedirs(path_to_write)
                #print(path_to_write + filename)
                write_to_netcdf(path_to_write + filename, latitude, longitude, altitude, frequency, density)
        i += 1


def process_file_kaliningrad(filename):
    base_dir = "/home/nick/code/python/netCDF"

    file_to_read = open(filename, 'r')
    lines = file_to_read.readlines()
    file_to_read.close()

    for line in lines:
        splitted_line = line.split()
        print(filename)
        date_str = '{}.{}.{}'.format(splitted_line[0], filename.split('/')[8], filename.split('/')[7])
        current_date = datetime.datetime.strptime(date_str, '%d.%b.%Y').date()
        print(current_date)

        path_to_write = '{}/kaliningrad/{}/{}.{:03d}/'.format(
            base_dir, current_date.year, current_date.year,
            current_date.timetuple().tm_yday
        )
        if not os.path.exists(path_to_write):
            os.makedirs(path_to_write)

        for hour in range(0, 24):
            if (
                splitted_line[hour + 1] == 'a' or
                splitted_line[hour + 1] == 'а' or
                splitted_line[hour + 1] == 'c' or
                splitted_line[hour + 1] == 'с' or
                splitted_line[hour + 1] == 'd' or
                splitted_line[hour + 1] == 'D' or
                splitted_line[hour + 1] == 'e' or
                splitted_line[hour + 1] == 'е' or
                splitted_line[hour + 1] == 'E' or
                splitted_line[hour + 1] == 'Е' or
                splitted_line[hour + 1] == 'b' or
                splitted_line[hour + 1] == 'B' or
                splitted_line[hour + 1] == 'В'
            ):
                continue

            current_time = datetime.datetime.strptime(
                '{:02d}0000'.format(hour), '%H%M%S'
            ).time()
            filename_to_save = make_filename(
                current_date.year, current_date.timetuple().tm_yday,
                current_time.hour, current_time.minute, 1, 1
            )

            latitude = 54
            longitude = 20
            altitude = 250
            frequency = float(splitted_line[hour + 1])
            density = (frequency / 0.124 / 100000) ** 0.5

            write_to_netcdf(
                path_to_write + filename_to_save, latitude, longitude,
                altitude, frequency, density
            )


def create_rozanov_file(filename):
    base_dir = "/home/nick/code/python/netCDF"

    file_to_read = open(filename, 'r')
    lines = file_to_read.readlines()
    file_to_read.close()

    for line in lines:
        splitted_line = line.split()
        print(filename)
        date_str = '{}.{}.{}'.format(splitted_line[0], filename.split('/')[8],
                                     filename.split('/')[7])
        current_date = datetime.datetime.strptime(date_str, '%d.%b.%Y').date()
        print(current_date)

        path_to_write = '{}/rozanov/'.format(base_dir)
        if not os.path.exists(path_to_write):
            os.makedirs(path_to_write)

        data_str = ''
        for hour in range(0, 24):
            if (
                splitted_line[hour + 1] == 'a' or
                splitted_line[hour + 1] == 'а' or
                splitted_line[hour + 1] == 'c' or
                splitted_line[hour + 1] == 'с' or
                splitted_line[hour + 1] == 'd' or
                splitted_line[hour + 1] == 'D' or
                splitted_line[hour + 1] == 'e' or
                splitted_line[hour + 1] == 'е' or
                splitted_line[hour + 1] == 'E' or
                splitted_line[hour + 1] == 'Е' or
                splitted_line[hour + 1] == 'b' or
                splitted_line[hour + 1] == 'B' or
                splitted_line[hour + 1] == 'В'
            ):
                continue

            time_hour = hour - 2
            if time_hour < 0:
                time_hour += 24
            current_time = '{}.00'.format(time_hour)
            f0f2 = float(splitted_line[hour + 1])
            data_str += '{} {} {}  {}  {}\n'.format(
                current_date.month, current_date.day,
                current_date.timetuple().tm_yday, current_time, f0f2
            )

        filename_to_save = make_rozanov_filename(
            current_date.year, current_date.timetuple().tm_yday
        )
        write_to_rozanov(path_to_write + filename_to_save, data_str)


def make_rozanov_filename(year, day_of_year):
    filename = 'MP{}{:03d}.DAT'.format(str(year)[-2:], day_of_year)

    return filename


def write_to_rozanov(filename, data):
    file_to_write = open(filename, 'a')

    headers_string = 'Mo da yda  time  FoF2\n'
    file_to_write.write(headers_string)
    file_to_write.write(data)

    file_to_write.close()


def fill_variables(file_netCDF, name, data, typ ='f4'):
    file_netCDF.createDimension(name, len(data))
    file_netCDF.createVariable(name, typ, (name, ))
    file_netCDF.variables[name][:] = data


def write_to_netcdf(filename, latitude, longitude, altitude, frequency, density):
    file_to_write = Dataset(filename, 'w', format='NETCDF4')
    f0f2 = (density ** 2) * 0.124 * 100000

    fill_variables(file_to_write, 'MSL_alt', [altitude, altitude, altitude])
    fill_variables(file_to_write, 'GEO_lat', [latitude, latitude, latitude])
    fill_variables(file_to_write, 'GEO_lon', [longitude, longitude, longitude])
    fill_variables(file_to_write, 'ELEC_dens', [f0f2, f0f2, f0f2])

    file_to_write.edmaxlat = latitude
    file_to_write.edmaxlon = longitude
    file_to_write.edmaxalt = altitude
    file_to_write.critfreq = frequency

    file_to_write.close()


def readFromNETCDF(filename, parameter_name):
    file_to_read = Dataset(filename, 'r', format = 'NETCDF3')
    # print(file_to_read)
    # MSL_alt = file_to_read.variables['MSL_alt'][:]
    # GEO_lat = file_to_read.variables['GEO_lat'][:]
    # GEO_lon = file_to_read.variables['GEO_lon'][:]
    # ELEC_dens = file_to_read.variables['ELEC_dens'][:]
    # edmaxlat = float(file_to_read.edmaxlat)
    # edmaxlon = float(file_to_read.edmaxlon)
    # edmaxalt = float(file_to_read.edmaxalt)
    # critfreq = float(file_to_read.critfreq)
    #
    # print(MSL_alt, GEO_lat, GEO_lon, ELEC_dens)
    # print(edmaxlat, edmaxlon, edmaxalt, critfreq)
    values = file_to_read.variables[parameter_name][:]
    file_to_read.close()

    return values


def make_filename(year, day, hour, minute, seans, number):
    file_name = 'ionPrf_C001.{}.{:03d}.{:02d}.{:02d}.{}.{}_nc'.format(year, day, hour, minute, seans, number)
    return file_name


def create_gsmtip_files_from_hammonia(file_to_read):
    start_date = datetime.strptime('2007-12-31 00:00:00', "%Y-%m-%d %H:%M:%S")
    times = readFromNETCDF(file_to_read, 'time')
    longitude_old_cells = readFromNETCDF(file_to_read, 'lon')
    latitude_old_cells = list(readFromNETCDF(file_to_read, 'lat'))
    latitude_old_cells.reverse()

    zonal_speed = readFromNETCDF(file_to_read, 'u')
    meridional_speed = readFromNETCDF(file_to_read, 'v')
    temperatures = readFromNETCDF(file_to_read, 'st')
    mm_air = readFromNETCDF(file_to_read, 'mm_air')

    for time_index in range(0, len(times)):
        time_delta = timedelta(days=times[time_index])
        corrected_date = start_date + time_delta
        day_of_year = corrected_date.timetuple().tm_yday
        current_hour = corrected_date.hour

        filename_zwind = 'GSMTIP/Zwind/Zwind_{}'.format(day_of_year)
        filename_vwind = 'GSMTIP/Vwind/Vwind_{}'.format(day_of_year)
        filename_temperature = 'GSMTIP/Temperature/T_{}'.format(day_of_year)
        filename_density = 'GSMTIP/Density/PL_{}'.format(day_of_year)

        file_to_write_zwind = open(filename_zwind, 'a')
        file_to_write_vwind = open(filename_vwind, 'a')
        file_to_write_temperature = open(filename_temperature, 'a')
        file_to_write_density = open(filename_density, 'a')
        file_to_write_zwind_surfer = open('{}_surfer'.format(filename_zwind), 'a')
        file_to_write_vwind_surfer = open('{}_surfer'.format(filename_vwind), 'a')
        file_to_write_temperature_surfer = open('{}_surfer'.format(filename_temperature), 'a')
        file_to_write_density_surfer = open('{}_surfer'.format(filename_density), 'a')

        surfer_header = 'longitude   latitude   value\n'
        file_to_write_zwind_surfer.write(surfer_header)
        file_to_write_vwind_surfer.write(surfer_header)
        file_to_write_temperature_surfer.write(surfer_header)
        file_to_write_density_surfer.write(surfer_header)

        for latitude in range(-85, 95, 5):
            latitude_i = int((90-latitude) / 5)
            header_message = '  UT={:10d} lat_i={:10d}\n   '.format(
                current_hour, latitude_i
            )

            file_to_write_temperature.write(header_message)
            file_to_write_density.write(header_message)
            file_to_write_zwind.write(header_message)
            file_to_write_vwind.write(header_message)

            latitude -= 2.5
            latitude_indexes = find_near_values(
                latitude_old_cells, latitude
            )

            for number, longitude in enumerate(range(5, 365, 5)):
                TMP_message = None
                longitude -= 2.5
                longitude_indexes = find_near_values(
                    longitude_old_cells, longitude
                )

                current_coordinates = [longitude, latitude]

                zonal_wind = bilinear_interpolation(
                    zonal_speed[time_index][0],
                    longitude_old_cells, latitude_old_cells,
                    longitude_indexes, latitude_indexes,
                    current_coordinates
                )
                meridional_wind = bilinear_interpolation(
                    meridional_speed[time_index][0],
                    longitude_old_cells, latitude_old_cells,
                    longitude_indexes, latitude_indexes,
                    current_coordinates
                )
                temperature = bilinear_interpolation(
                    temperatures[time_index][0],
                    longitude_old_cells, latitude_old_cells,
                    longitude_indexes, latitude_indexes,
                    current_coordinates
                )
                TMP = bilinear_interpolation(
                    mm_air[time_index][0],
                    longitude_old_cells, latitude_old_cells,
                    longitude_indexes, latitude_indexes,
                    current_coordinates
                )
                density = calculate_air_density(TMP, temperature)

                temperature_data = '  {:.4E}'.format(temperature)
                density_data = '  {:.4E}'.format(density)
                zwind_data = '  {:.4E}'.format(zonal_wind)
                vwind_data = '  {:.4E}'.format(-meridional_wind)

                longitude_i = number + 1
                if longitude_i % 6 == 0:
                    TMP_message = '\n'
                    if longitude_i != 72:
                        TMP_message += '   '

                file_to_write_temperature.write(temperature_data)
                file_to_write_density.write(density_data)
                file_to_write_zwind.write(zwind_data)
                file_to_write_vwind.write(vwind_data)

                if TMP_message:
                    file_to_write_temperature.write(TMP_message)
                    file_to_write_density.write(TMP_message)
                    file_to_write_zwind.write(TMP_message)
                    file_to_write_vwind.write(TMP_message)

                surfer_coordinates = '\n{}   {} '.format(longitude, latitude)
                file_to_write_zwind_surfer.write(surfer_coordinates)
                file_to_write_vwind_surfer.write(surfer_coordinates)
                file_to_write_temperature_surfer.write(surfer_coordinates)
                file_to_write_density_surfer.write(surfer_coordinates)

                file_to_write_zwind_surfer.write(zwind_data)
                file_to_write_vwind_surfer.write(vwind_data)
                file_to_write_temperature_surfer.write(temperature_data)
                file_to_write_density_surfer.write(density_data)


def calculate_air_density(mm_air, temperature):
    density = mm_air * PRESSURE / temperature / GAS_CONSTANT / 1000000
    # print(mm_air, temperature, density)

    return density


def bilinear_interpolation(
    values, longitude_cells, latitude_cells,
    longitude_indexes, latitude_indexes, current_coordinates
):

    longitude = current_coordinates[0]
    latitude = current_coordinates[1]
    top_left = values[latitude_indexes[1]][longitude_indexes[0]]
    top_right = values[latitude_indexes[1]][longitude_indexes[1]]
    bottom_left = values[latitude_indexes[0]][longitude_indexes[0]]
    bottom_right = values[latitude_indexes[0]][longitude_indexes[1]]

    denominator = longitude_cells[longitude_indexes[1]] - longitude_cells[longitude_indexes[0]]
    numerator_right_side = longitude_cells[longitude_indexes[1]] - longitude
    numerator_left_side = longitude - longitude_cells[longitude_indexes[0]]

    interpolated_top = numerator_right_side / denominator * top_left
    interpolated_top += numerator_left_side / denominator * top_right
    interpolated_bottom = numerator_right_side / denominator * bottom_left
    interpolated_bottom += numerator_left_side / denominator * bottom_right

    denominator = latitude_cells[latitude_indexes[1]] - latitude_cells[latitude_indexes[0]]
    numerator_top_side = latitude_cells[latitude_indexes[1]] - latitude
    numerator_bottom_side = latitude - latitude_cells[latitude_indexes[0]]

    converted_value = numerator_top_side / denominator * interpolated_bottom
    converted_value += numerator_bottom_side / denominator * interpolated_top

    return converted_value


def find_near_values(old_cells, coordinate):
    for index in range(0, len(old_cells) - 1):
        if old_cells[index] <= coordinate <= old_cells[index + 1]:
            return index, index + 1

    return 0, 0


# all_month_names = [
#     'jan', 'feb', 'mar', 'apr', 'may', 'jun',
#     'jul', 'aug', 'sep', 'oct', 'nov', 'dec'
# ]
#
# print(makeFileName(2007, 1, 3, 5, 1234, 123))
# processFile('Winter.dat')
# for year in range(2007, 2017):
#     base_path = '/home/nick/code/python/netCDF/convert_data/{}/'.format(year)
#     for month_name in all_month_names:
#         path = base_path + month_name
#         # process_file_kaliningrad(path)
#         create_rozanov_file(path)

# path_to_hammonia = '/media/nick/Elements/PanOply/ham_test2_200901.nc'
path_to_hammonia = '/media/nick/Elements/PanOply/HAMMONIA_test_200901_6h.nc'
create_gsmtip_files_from_hammonia(path_to_hammonia)

# processFile('Summer.dat')
# writeToNETCDF('temp')
# readFromNETCDF('ham_test2_200901.nc')
