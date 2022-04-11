import pandas as pd
import argparse
import sys
import time
import os
import shutil


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
    return ("\nHow to run? It's easy, just follow the terminal examples listed.\n"
            "\nTo run feature extraction only, run the following commands:"
            "\npython MALDI_TOF_features.py -f directory_with_raw_maldi_csv_files \n"
            "\nTo run feature extraction and antibiotic label matching, run the following:"
            "\npython MALDI_TOF_features.py --folder directory_with_raw_maldi_csv_files --labels labels_dir\n"
)


parser = argparse.ArgumentParser(description='MALDI-TOF MS 3 dalton feature extraction binning from bacterial spectra samples'
                                             '\nFeatures are extracted from the 2-20 kdaltons m/z range',
                                 epilog=_epilog(),
                                 formatter_class=argparse.RawDescriptionHelpFormatter)

parser.add_argument(
    "-f",
    "--folder",
    help="MALDI spectra folder",
    type=str,
    required=True)

parser.add_argument(
    "-o",
    "--feature_output_dir",
    help="Folder output name",
    type=str,
    default='MALDI-TOF-MS Features')

parser.add_argument(
    '-l',
    '--labels',
    help='labels output folder',
    type=str
)


args = parser.parse_args()

# navigate into the user specified directory
original_path = os.getcwd()


def extraction():
    typewrite('Extracting features from MALDI-TOF MS spectra\n')
    os.chdir(args.folder)
    spectra_files = os.listdir(".")
    spectra_files = sorted(spectra_files)
    dataframe_dict = {}

    s = list(range(1000, 10000, 3))
    e = list(range(1003, 10003, 3))

    names = []
    for x, y in zip(s, e):
        name = '{}-{}'.format(str(x), str(y))
        names.append(name)

    for spectrum in spectra_files:
        # condition for spectra file, other file extensions ignored
        if spectrum.endswith('.csv'):
            typewrite('Processing {}... \n'.format(spectrum.split('.')[0]))

            df = pd.read_csv(spectrum)

            if not df.shape[1] == 2:
                raise MaldiTofMsException(df.shape[1], spectrum)

            # example spectra which has been trimmed to 10_000 m/z

            stop_value = 9999.734225
            start_spectra_val = 1000
            end_spectra_val = 1003
            intensity_m_z_vals = []

            for i in range(len(df)):

                if start_spectra_val < stop_value:

                    df_spectra = df[df['mass'].between(start_spectra_val, end_spectra_val, inclusive='left')]
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

        else:
            continue

    # creating a directory for MALDI-TOF-MS feature output
    try:
        os.mkdir(args.feature_output_dir)
    except FileExistsError:
        pass


    dataframe.to_csv("MALDI_features.csv")
    shutil.move("MALDI_features.csv", args.feature_output_dir)

    # if the folder already exists, remove it and create another

    maldi_dir = original_path + '/' + args.feature_output_dir
    os.chdir(original_path)
    if os.path.exists(maldi_dir):
        shutil.rmtree(maldi_dir)
    os.chdir(args.folder)
    shutil.move(args.feature_output_dir, original_path)


    typewrite("Complete...\nMALDI features now in a csv file called MALDI_features.csv "
              "in a directory called {}\n".format(args.feature_output_dir))


# part 2 of script
# navigate to the folder with the antibiotic label information

def lables_match():
    typewrite('\nMatching antibiotic resistance labels with features\n')
    labels_dir = original_path + '/' + str(args.labels)
    os.chdir(labels_dir)

    features_csv = 'MALDI_features.csv'
    features_df = pd.read_csv(original_path + '/' + str(args.feature_output_dir) + '/' + features_csv, index_col=0)

    antibiotic_labels_dict = {}
    for label in sorted(os.listdir('.')):
        if label.endswith('csv'):
            df = pd.read_csv(label)
            name = label.split('.')[0]
            df.index = [name]
            antibiotic_labels_dict[name] = df

    antibiotic_df = pd.concat([i for i in antibiotic_labels_dict.values()])
    feat_and_labels = features_df.join(antibiotic_df)
    feat_and_labels.to_csv('MALDI_TOF_features_and_labels.csv')

    if os.path.exists('Features_and_labels_ML_dir'):
        shutil.rmtree('Features_and_labels_ML_dir')
    os.mkdir('Features_and_labels_ML_dir')

    if os.path.exists(original_path + '/' + 'Features_and_labels_ML_dir'):
        shutil.rmtree(original_path + '/' + 'Features_and_labels_ML_dir')

    shutil.move('MALDI_TOF_features_and_labels.csv', 'Features_and_labels_ML_dir')
    shutil.move(os.getcwd() + "/" + 'Features_and_labels_ML_dir', original_path)

    typewrite('\nFeatures matched with antibiotic labels complete\n')

    typewrite("\nPlease cite https://github.com/StephenFordham 2022, if you use this script\n")


# calling logic
if args.folder and not args.labels:
    extraction()
elif args.folder and args.labels:
    extraction()
    lables_match()





