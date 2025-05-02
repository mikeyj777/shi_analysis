from data.db_connection import get_data_table_return_as_dataframe

# get all rows from "HazardStudy" db table - index = StudyID

# haz_studies_df = get_data_table_return_as_dataframe("HazardStudy")

# get all rows from "StudyChemicalComposition" table - index = HazardStudy
# chem_comps_df = get_data_table_return_as_dataframe("StudyChemicalComposition")

# parse inputs to run each tt

#simulating for now until back on emn connection
from py_lopa.model_interface import Model_Interface

m_io = Model_Interface

# store shi at process and storage conditions (based on key comp for toxic and overall for flammable)

# store footprints for lfl / 50% lfl / 25% lfl / erpg-1 / erpg-2 / erpg-3 / 10xerpg-3