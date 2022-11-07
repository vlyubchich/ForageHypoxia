import pandas as pd
biomass = r'C:\Users\Manke\Documents\GitHub\ForageHypoxia\data_benthos\benthos_biomass_RANDOM.csv'
taxa = r'C:\Users\Manke\Documents\GitHub\ForageHypoxia\data_benthos\Taxa_group_IDs_updated_rjw_v2.xlsx - Invert_IDs.csv'
df = pd.read_csv(biomass, encoding = 'ANSI')
trjw = pd.read_csv(taxa, encoding = 'ANSI').astype(str)
trjw = trjw.astype(str) #this is here because the order column is classified as integers as opposed to strings

'''
Taxa sorting functions
'''
def count_list(sorted_list, set_list):
    the_list = []
    for i in set_list:
        the_list.append(sorted_list.count(i))
    return the_list

def taxa_spreadsheet(dataframe, column):
    data = {column: sorted(list(set(dataframe[column]))),
            'Count': count_list(sorted(list(dataframe[column])), sorted(list(set(dataframe[column]))))}
    spreadsheet = pd.DataFrame(data) #makes the spreadsheet
    return spreadsheet

#creates spreadsheet
def make_spreadsheet(sheet, name):
    file_name = name + '.csv'
    sheet.to_csv(file_name)

'''
Compiling Taxa by function groups, aggregate groups, class, order, and family
How many times they occur in a sample
'''
taxdat = taxa_spreadsheet(df, 'LBL') #taxa data
taxlist = list(trjw.columns) #lists the different grouping columns (family, order, etc.)
for col in taxlist: #col is short for column
    taxa_spreadsheet(trjw, col)

print(taxa_spreadsheet(trjw, 'Family'))