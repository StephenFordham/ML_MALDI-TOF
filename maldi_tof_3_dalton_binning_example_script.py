import pandas as pd


df = pd.read_csv('spectrum1.csv')

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


