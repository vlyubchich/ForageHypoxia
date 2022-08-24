import pandas as pd
biomass = r'C:\Users\Manke\Documents\GitHub\ForageHypoxia\data_benthos\benthos_biomass_RANDOM.csv'
df = pd.read_csv(biomass, encoding = 'ANSI')

'''
Creates spreadsheet
'''
def make_spreadsheet(sheet, name):
    file_name = name + '.csv'
    sheet.to_csv(file_name)
'''
Functions for retrieving sum, mean, and count
var = variable
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
creates spreadsheet
'''
class spread:
    def __init__(self, var):
        self.sum = retrieve_biomass_sum(var)
        self.mean = retrieve_biomass_mean(var)
        self.count = retrieve_biomass_count(var)

    def final_spread(self):
        return pd.concat([self.sum, self.mean, self.count], axis=1)

annual_data = spread('Year').final_spread()

'''
Counts number of any given variable in a year, in this instance, number of stations per year
var1 = 
var2 = 
y = 
'''
def variable_count(var1, var2, y):
    simp_df = df[df[var1] == y] #simp_df stands for simplified dataframe
    return len(set(simp_df[var2]))

def variable_count_column(var1, var2):
    count_list = []
    for y in list(set(df[var1])):
        count_list.append(variable_count(var1, var2, y))
    return count_list

#Couldn't figure out how to rename COUNT column, so instead I just made a new column with the same values, and deleted COUNT
#Species_Count is the number of organisms that were gathered that year
annual_data['Species_Count'] = annual_data['COUNT']
annual_data.drop('COUNT', axis = 1, inplace = True)
#Species_Total is the number of species identified that year (not total number of organisms
annual_data['Species_Total'] = variable_count_column('Year', 'LBL')
#Station_Count is the number of stations used that year
annual_data['Station_Count'] = variable_count_column('Year', 'STATION')
print(annual_data)

#make_spreadsheet(annual_data, 'Annual_Data')