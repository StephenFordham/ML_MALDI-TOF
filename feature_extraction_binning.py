import pandas as pd
import argparse
import sys
import time


# terminal output retro, slow function
def typewrite(message):
    for char in message:
        sys.stdout.write(char)
        sys.stdout.flush()

        if char != '\n':
            time.sleep(0.05)
        else:
            time.sleep(0.5)


# spectra file input validation
def spectra_input_validation(spectra_csv_file):
    if not spectra_csv_file.endswith('.csv'):
        raise argparse.ArgumentTypeError('Please enter a valid csv file')
    df = pd.read_csv(spectra_csv_file)
    if df.shape[1] != 2:
        raise argparse.ArgumentTypeError('Please enter a spectra file with only m/z and intensity values')
    return spectra_csv_file


parser = argparse.ArgumentParser()

parser.add_argument(
                    "-s",
                    "--spectra_file",
                    help="Please enter the path to your spectra csv file",
                    type=spectra_input_validation,
                    required=True)


parser.add_argument(
                    "-f",
                    "--filename",
                    help="File output name",
                    type=str,
                    default="spectra_file")


args = parser.parse_args()

typewrite('Extracting features from MALDI-TOF MS spectra\n')


df = pd.read_csv(args.spectra_file)
s = list(range(1000, 10000, 3))
e = list(range(1003, 10003, 3))

names = []
for x, y in zip(s,e):
    name = '{}-{}'.format(str(x), str(y))
    names.append(name)

# example spectra which has been trimmed to 10_000 m/z

stop_value = 9999.734225
start_spectra_val = 1000
end_spectra_val = 1003
intensity_m_z_vals = []


for i in range(len(df)):

    if not start_spectra_val > stop_value:

        df_spectra = df[df['mass'].between(start_spectra_val, end_spectra_val)]
        intensity_3_dalton_window_vals = df[df['mass'].between(start_spectra_val, end_spectra_val)]['intensity'].sum()
        intensity_m_z_vals.append(intensity_3_dalton_window_vals)
        start_spectra_val += 3
        end_spectra_val += 3

    else:
        break

# The dataframe represents a single MALDI-TOF sample

maldi_tof_df = pd.DataFrame(data=intensity_m_z_vals).T
maldi_tof_df.columns = [names]

maldi_tof_df.to_csv(args.filename + ".csv")

typewrite("Complete...\nOutput file in directory\n")