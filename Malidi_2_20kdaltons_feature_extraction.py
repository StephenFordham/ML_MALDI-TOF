import pandas as pd
import argparse
import sys
import time
import os
import shutil
import subprocess

class MaldiTofMsException(Exception):
    def __init__(self, number_of_cols, spectrum):
        self._col_num = number_of_cols
        self._spectrum = spectrum

    def __str__(self):
        return '\n{} columns detected for {}\n' \
               'Only M/Z and intensity values are valid columns'.format(self._col_num, self._spectrum)


# terminal output retro, slow function
def typewrite(message):
    for char in message:
        sys.stdout.write(char)
        sys.stdout.flush()

        if char != '\n':
            time.sleep(0.05)
        else:
            time.sleep(0.5)

def _epilog():
    return ("\n How to run? It's easy, just follow the terminal examples listed"
            "\n python Malidi_2_20kdaltons_feature_extraction.py --folder directory_with_raw_maldi_csv"
            "\n python Malidi_2_20kdaltons_feature_extraction.py -f directory_with_raw_maldi_csv")


parser = argparse.ArgumentParser(description='MALDI-TOF MS 3 dalton feature extraction from bacterial spectra samples'
                                             '\nFeatures are extracted from 2-20 kdaltons',
                                 epilog=_epilog(),
                                 formatter_class=argparse.RawDescriptionHelpFormatter)

parser.add_argument(
    "-f",
    "--folder",
    help="File output name",
    type=str)

parser.add_argument(
    "-o",
    "--output_folder",
    help="Folder output name",
    type=str,
    default='MALDI-TOF-MS Features')


args = parser.parse_args()

typewrite('Extracting features from MALDI-TOF MS spectra\n')

# navigate into the user specified directory
original_path = os.getcwd()
os.chdir(args.folder)
spectra_files = os.listdir('.')
dataframe_dict = {}

for spectrum in spectra_files:
    print(spectrum)

    df = pd.read_csv(spectrum)

    if not df.shape[1] == 2:
        raise MaldiTofMsException(df.shape[1], spectrum)

    s = list(range(1000, 10000, 3))
    e = list(range(1003, 10003, 3))

    names = []
    for x, y in zip(s, e):
        name = '{}-{}'.format(str(x), str(y))
        names.append(name)

    # example spectra which has been trimmed to 10_000 m/z

    stop_value = 9999.734225
    start_spectra_val = 1000
    end_spectra_val = 1003
    intensity_m_z_vals = []

    for i in range(len(df)):

        if start_spectra_val < stop_value:

            df_spectra = df[df['mass'].between(start_spectra_val, end_spectra_val)]
            intensity_3_dalton_window_vals = df_spectra['intensity'].sum()
            intensity_m_z_vals.append(intensity_3_dalton_window_vals)
            start_spectra_val += 3
            end_spectra_val += 3

        else:
            break

    # The dataframe represents a single MALDI-TOF sample

    maldi_tof_df = pd.DataFrame(data=intensity_m_z_vals).T
    maldi_tof_df.columns = [names]
    maldi_tof_df.index = [spectrum.split('.')[0]]

    # update dictionary with spectrum name and spectrum dataframe on each iteration of the for loop
    dataframe_dict[str(spectrum)] = maldi_tof_df
    dataframe = pd.concat([value for value in dataframe_dict.values()])

# creating a directory for MALDI-TOF-MS feature output
directory = args.output_folder
os.mkdir(directory)


dataframe.to_csv("MALDI_features.csv")
shutil.move("MALDI_features.csv", directory)
shutil.move(os.getcwd() + "/" + directory, original_path)


typewrite("Complete...\nMALDI features now in a csv file called MALDI_features.csv "
          "in a directory called {}\n".format(args.output_folder))
typewrite("\nPlease cite https://github.com/StephenFordham 2022, if you use this script\n")