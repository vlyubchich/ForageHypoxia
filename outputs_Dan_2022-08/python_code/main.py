import pandas as pd
biomass = r'C:\Users\Manke\Documents\GitHub\ForageHypoxia\data_benthos\benthos_biomass_RANDOM.csv'
taxa = r'C:\Users\Manke\Documents\GitHub\ForageHypoxia\data_benthos\Taxa_group_IDs_updated_rjw_v2.xlsx - Invert_IDs.csv'
df = pd.read_csv(biomass, encoding = 'ANSI') #df stands for dataframe
df['SAMPLE_DATE'] = df['SAMPLE_DATE'].apply(pd.to_datetime)
trjw = pd.read_csv (taxa, encoding = 'ANSI')

'''
Functions for biomass by station, given separate functions because I needed other variables to be included in the spreadsheet, as seen below (especially coordinates).
However, I recently learned how to add the columns for those variables without needing to make entirely separate functions, have not made that change due to time constraints.
'''
#calculates total biomass
def retrieve_station_biomass_sum(var):
    grouped_df = df.groupby([var, 'CBPSEG', 'cbpseg_st', 'LATITUDE', 'LONGITUDE']).sum()
    grouped_df['SUM'] = grouped_df['VALUE'] #renames the column
    return grouped_df['SUM']

#calculates mean biomass
def retrieve_avg_station_biomass(var):
    grouped_df = df.groupby([var, 'CBPSEG', 'cbpseg_st', 'LATITUDE', 'LONGITUDE']).mean()
    grouped_df['MEAN'] = grouped_df['VALUE']
    return grouped_df['MEAN']

def retrieve_station_biomass_count(var):
    grouped_df = df.groupby([var, 'CBPSEG', 'cbpseg_st', 'LATITUDE', 'LONGITUDE']).count()
    grouped_df['COUNT'] = grouped_df['VALUE']
    return grouped_df['COUNT']

'''
This function creates the spreadsheet
'''
def make_spreadsheet(sheet, name):
    file_name = name + '.csv'
    sheet.to_csv(file_name)

'''
Biomass for CBPSEG, cbpseg_site, and other functions. These are the standard functions you will use for every other column besides station.
'''
def retrieve_biomass_sum(var):
    the_sum = df.groupby([var]).sum()
    the_sum['SUM'] = the_sum['VALUE']  # renames the column
    return the_sum['SUM']
def retrieve_biomass_mean(var):
    the_mean = df.groupby([var]).mean()
    the_mean['MEAN'] = the_mean['VALUE']  # renames the column
    return the_mean['MEAN']
def retrieve_biomass_count(var):
    the_count = df.groupby([var]).count()
    the_count['COUNT'] = the_count['VALUE']  # renames the column
    return the_count['COUNT']

'''
Creates spreadsheets
'''
class station_spread: #specifically for Stations
    def __init__(self, var):
        self.sum = retrieve_station_biomass_sum(var)
        self.mean = retrieve_avg_station_biomass(var)
        self.count = retrieve_station_biomass_count(var)

    def final_spread(self):
        return pd.concat([self.sum, self.mean, self.count], axis=1)

class spread:
    def __init__(self, var):
        self.sum = retrieve_biomass_sum(var)
        self.mean = retrieve_biomass_mean(var)
        self.count = retrieve_biomass_count(var)

    def final_spread(self):
        return pd.concat([self.sum, self.mean, self.count], axis=1)


Station_Data = station_spread('STATION').final_spread()
CBPSEG_Data = spread('CBPSEG').final_spread()

print(Station_Data)

'''
make_spreadsheet(Station_Data, 'Station_Data')
make_spreadsheet(CBPSEG_Data, 'CBPSEG_Data')
make_spreadsheet(cbpseg_st_Data, 'cbpseg_st_Data')
make_spreadsheet(date_Data, 'date_Data')
make_spreadsheet(species_Data, 'species_Data')
'''

print('Spreadsheets complete!')