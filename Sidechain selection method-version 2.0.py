#!/usr/bin/env python
# coding: utf-8

# # Part1: construct a receptor-peptide residue interaction library

# In[6]:


import pandas as pd
import numpy as np
import requests

# retrieve all available structures 
url = "https://gpcrdb.org/services/structure/"
response = requests.get(url)
data_strucutre = response.json()
df_all_structure = pd.DataFrame(data_strucutre)

# identifier for protein structure
# class == class A & ligands == peptides or proteins
class_category = df_all_structure['class'].value_counts() # class_category
mask = df_all_structure.loc[:,'class'].str.startswith('Class A') 
df_ligands_type = df_all_structure['ligands'].apply(lambda x: x[0]['type'] if x else None) # extract the ligand type
ligands_type_category = df_ligands_type.value_counts()
mask_2 = df_ligands_type.isin(["peptide", "protein"])
df_all_structure = df_all_structure[mask & mask_2]

# remove repeated *protein - ligands_name* pair and then keep the best resolution
df_ligands_name = df_all_structure['ligands'].apply(lambda x: x[0]['name'] if x else None) # extract the ligand name
df_ligands_name = df_ligands_name.rename('ligands_name') 

# # add ligand name column into the dataframe
df_full_structure_with_ligands_name = pd.concat([df_all_structure, df_ligands_name], axis=1)
best_resolution_index = df_full_structure_with_ligands_name.groupby(['protein', 'ligands_name'])['resolution'].idxmin()
filter_with_identifier = df_all_structure.loc[best_resolution_index]

# get the pdb_code that fulfil the condition
df_structure = filter_with_identifier
pdb_code_list = df_structure["pdb_code"].tolist()

# get all the interaction data for each structure complex
all_interaction_data = []

for code in pdb_code_list:
    url = f"https://jimmy.gpcrdb.org/services/structure/{code}/peptideinteraction/" # Get a list of interactions between structure and peptide ligand
    response = requests.get(url)
    if response.status_code == 200:
        all_interaction_data.extend(response.json())
    else:
        print(f"Failed to fetch data from {url}")
raw_interaction =  pd.DataFrame(all_interaction_data)

# keep those equal to sidechain sidechain interaction/Sidenchain Backbone groups
df_SS_SB_interaction = raw_interaction[(raw_interaction['structural_interaction'] == 'SS') | (raw_interaction['structural_interaction'] == 'SB')]

# remove the wan-der-wons interaction
df_SS_SB_interaction_without_van_der_waals = df_SS_SB_interaction[df_SS_SB_interaction["interaction_type"] != "van-der-waals"]

# removing peptide_amnio_acid =='X'
df_SS_SB_interaction_without_van_der_waals = df_SS_SB_interaction_without_van_der_waals[df_SS_SB_interaction_without_van_der_waals["peptide_amino_acid"] != "X"]

# ignore the interaction type ANDT keep the unique pdb_code-Generic number-peptide AA-receptor AA pair # different pdb_code
unique_structure_df_interaction = df_SS_SB_interaction_without_van_der_waals.drop_duplicates(subset=['pdb_code', 'receptor_residue_generic_number', 'receptor_amino_acid', 'peptide_amino_acid'])

# filter 1
df_interaction = unique_structure_df_interaction.loc[:, ["receptor_residue_generic_number",  "receptor_amino_acid", "peptide_amino_acid", "ca_distance", "ca_cb_angle"]]

# convert the ca_cb_angle and ca_ca_distance column to float unless return NaN                                              
df_interaction["ca_cb_angle"] = pd.to_numeric(df_interaction["ca_cb_angle"], errors='coerce') 
df_interaction["ca_distance"] = pd.to_numeric(df_interaction["ca_distance"], errors='coerce')      

# find the rows with NaN value
rows_with_nan = df_interaction[df_interaction.isna().any(axis=1)]

# remove rows with NaN value(must have ca_distance and ca_cb_angle)
df_interaction = df_interaction.dropna(subset=['ca_distance', 'ca_cb_angle'])
 
# sort in order
df_interaction = df_interaction.sort_values(["receptor_residue_generic_number","receptor_amino_acid"])

# obtain median and sd value for each RGN - receptor - AA unqiue pair 
df1 = df_interaction.groupby(['receptor_residue_generic_number', 'receptor_amino_acid', 'peptide_amino_acid']).agg(
    Distance_median =('ca_distance', 'median'),
    Distance_sd =('ca_distance', 'std'),
    Angle_median =('ca_cb_angle', 'median'),
    Angle_sd =('ca_cb_angle', 'std')
).reset_index()


# # obtain the count column for each RGN - receptor - AA unqiue pair 
df2 = df_interaction.groupby(['receptor_residue_generic_number', 'receptor_amino_acid', 'peptide_amino_acid']).size().reset_index(name='count')

# # combine together
df_merged = pd.merge(df1, df2,
                     on=['receptor_residue_generic_number', 'receptor_amino_acid', 'peptide_amino_acid'],
                     how='left')

# # obtain the sum column based on the count column 
df_merged['sum_total_count'] = df_merged.groupby(['receptor_residue_generic_number', 'receptor_amino_acid'])['count'].transform('sum')
df_merged['frequency(%)'] = (df_merged['count'] / df_merged['sum_total_count']) * 100

# Restrict decimals
# Rounding 'receptor_residue_generic_number' to two decimal places
df_merged['Distance_median'] = df_merged['Distance_median'].round(1)
df_merged['Distance_sd'] = df_merged['Distance_sd'].round(1)
df_merged['Angle_median'] = df_merged['Angle_median'].round(1)
df_merged['Angle_sd'] = df_merged['Angle_sd'].round(1)

# Rounding 'frequency(%)' to the nearest whole number and converting to integers
df_merged['frequency(%)'] = df_merged['frequency(%)'].round().astype(int)

reference_GPCR= df_merged.reset_index(drop=True)

reference_GPCR


# # Part2: Validation

# In[8]:


# 2.1 Require functions

# amino acids dictionary in Upper letter
amino_acids_dict = {
        "ALA": "A", "ARG": "R", "ASN": "N", "ASP": "D", "CYS": "C", "GLN": "Q", "GLU": "E",
        "GLY": "G", "HIS": "H", "ILE": "I", "LEU": "L", "LYS": "K",
        "MET": "M", "PHE": "F", "PRO": "P", "SER": "S",
        "THR": "T", "TRP": "W", "TYR": "Y", "VAL": "V",
        "X":"X"
    }

# vector calculation
def calculate_vector(row): 
    # v1:v2 = CA : CB = CA→CB vector
    v1 = np.array([row['x_coord_x'], row['y_coord_x'], row['z_coord_x']])
    v2 = np.array([row['x_coord_y'], row['y_coord_y'], row['z_coord_y']])

    # check nan
    if np.isnan(v1).any() or np.isnan(v2).any():
        return np.nan
    
    # v1 → v2
    V = v2 - v1

    # normalized to unit direction vector
    norm_V = np.linalg.norm(V)
    if norm_V != 0:
        V = V / norm_V
    else:
        V = None

    return V

# angle calculation
def calculate_angle(row): # based on two vector
    # extract the vector
    v1 = np.array(row["receptor_vector_angle"]) 
    v2 = np.array(row["ligand_vector_angle"])
    
    # calculate the dot product
    dot_product = np.dot(v1, v2)

    # calculate the magnitude
    norm_v1 = np.linalg.norm(v1)
    norm_v2 = np.linalg.norm(v2)

    # avoid dividing zero
    if norm_v1 == 0 or norm_v2 == 0:
        return None

    # Calculate the angle in radians
    cos_radian = dot_product / (norm_v1 * norm_v2)
    # limit the radian from 0 to π
    radian = np.arccos(np.clip(cos_radian, -1.0 , 1.0)) 
    # transform to angle degree
    angle = np.degrees(radian)

    return angle

# pre-treatment for the input data
def pre_treatment(inference_df, df_protein_residue):
    # combine the GRN with the inference_df
    merged = pd.merge(inference_df, df_protein_residue,
                     on = "residue_number",
                     how = "left" )

    # filter out those atom_name == CA/CB and GRN are not None
    mask_3 = (merged.loc[:, "atom_name"].isin(["CA", "CB"]) & (merged.loc[:, 'display_generic_number'].notnull()))
    merged_df = merged[mask_3]
    merged_df = merged_df.sort_values(by='residue_number', ascending=True)

    ligand_atoms = inference_df[inference_df["chain_id"] == 'B'].copy()
    ligand_atoms = ligand_atoms.loc[ligand_atoms["atom_name"].isin(["CA", "CB"])]
    # all residue turn to Ala
    # ligand_atoms.drop(["residue_name"], axis = 1 ,inplace= True)

    receptor_atoms = merged_df[merged_df["chain_id"].isin(['A', 'R'])] 
    
    # Check if receptor DataFrame is empty
    if receptor_atoms.empty:
        raise ValueError("Error warning! NO online GRN MATCH for this receptor")
    
    # filter the CA atom and convert to single letter
    receptor_atoms_CA = receptor_atoms[receptor_atoms.loc[:,"atom_name"] == "CA"].copy()
    ligand_atoms_CA = ligand_atoms[ligand_atoms.loc[:,"atom_name"] == "CA"].copy()

    # create the single_letter column 
    receptor_atoms_CA.loc[:,'receptor_amino_acid'] = receptor_atoms['residue_name'].map(amino_acids_dict)
    receptor_atoms_CA.drop(["residue_name","chain_id","residue_number","atom_name"], axis = 1, inplace = True)

    # keep the residue number for ligand
    ligand_atoms_CA.drop(["chain_id","atom_name"], axis = 1, inplace = True)

    return ligand_atoms ,receptor_atoms ,receptor_atoms_CA , ligand_atoms_CA, merged

# calculate the Ca-Ca dist
def calculate_CaCa_distance(receptor_atoms_CA, ligand_atoms_CA):
    # CA-CA distance
    # distance = np.sqrt((x1 - x0)**2 + (y1 - y0)**2 + (z1 - z0)**2)
    receptor_atoms_CA['key'] = 1
    ligand_atoms_CA['key'] = 1
    cross_joined_merged = pd.merge(receptor_atoms_CA, ligand_atoms_CA, on='key').drop('key', axis=1)

    cross_joined_merged['Distance'] = np.sqrt(
        (cross_joined_merged['x_coord_x'] - cross_joined_merged['x_coord_y']) ** 2 +
        (cross_joined_merged['y_coord_x'] - cross_joined_merged['y_coord_y']) ** 2 +
        (cross_joined_merged['z_coord_x'] - cross_joined_merged['z_coord_y']) ** 2
    )
    cross_joined_merged = cross_joined_merged.rename(columns={"display_generic_number": "receptor_residue_generic_number"})

    Ca_Ca_distance = cross_joined_merged[["receptor_residue_generic_number", "residue_number","receptor_amino_acid","Distance"]]
    return Ca_Ca_distance

# calculate the nominally existing Cb atom coordinates for GLY residue
from math import *
class Vector(tuple):
    """Tuple subclass implementing basic 3D vectors"""

    def __new__(cls, x, y, z):
        return tuple.__new__(cls, (float(x), float(y), float(z)))

    def perp(self, other):
        """Part perpendicular to other vector"""
        dp = self[0]*other[0] + self[1]*other[1] + self[2]*other[2]
        return Vector(self[0] - dp*other[0],
                      self[1] - dp*other[1],
                      self[2] - dp*other[2])

    @property
    def unit(self):
        """Scaled to unit length"""
        n = sqrt(self[0]*self[0] + self[1]*self[1] + self[2]*self[2])
        return Vector(self[0]/n, self[1]/n, self[2]/n)

    @property
    def norm(self):
        """Euclidean length"""
        return sqrt(self[0]*self[0] + self[1]*self[1] + self[2]*self[2])

    @property
    def normsqr(self):
        """Euclidean length squared"""
        return self[0]*self[0] + self[1]*self[1] + self[2]*self[2]

    @property
    def x(self):
        """Vector x coordinate"""
        return self[0]

    @property
    def y(self):
        """Vector y coordinate"""
        return self[1]

    @property
    def z(self):
        """Vector z coordinate"""
        return self[2]

    def __bool__(self):
        """Nonzero vector"""
        return (self[0]*self[0] + self[1]*self[1] + self[2]*self[2] > 0)

    def __abs__(self):
        """abs(a): Euclidean length of vector a"""
        return sqrt(self[0]*self[0] + self[1]*self[1] + self[2]*self[2])

    def __add__(self, other):
        """a + b: Vector addition"""
        if isinstance(other, (tuple, list)) and len(other) >= 3:
            return Vector(self[0]+other[0], self[1]+other[1], self[2]+other[2])
        else:
            return NotImplemented

    def __radd__(self, other):
        """b + a: Vector addition"""
        if isinstance(other, (tuple, list)) and len(other) >= 3:
            return Vector(other[0]+self[0], other[1]+self[1], other[2]+self[2])
        else:
            return NotImplemented

    def __mul__(self, other):
        """a * b: Scalar multiplication"""
        if isinstance(other, (int, float)):
            return Vector(self[0]*other, self[1]*other, self[2]*other)
        else:
            return NotImplemented

    def __rmul__(self, other):
        """b * a: Scalar multiplication"""
        if isinstance(other, (int, float)):
            return Vector(other*self[0], other*self[1], other*self[2])
        else:
            return NotImplemented

    def __neg__(self):
        """-a: Negation"""
        return Vector(-self[0], -self[1], -self[2])

    def __or__(self, other):
        """a | b: Dot product"""
        if isinstance(other, (tuple, list)) and len(other) >= 3:
            return self[0]*other[0] + self[1]*other[1] + self[2]*other[2]
        else:
            return NotImplemented

    def __ror__(self, other):
        """b | a: Dot product"""
        if isinstance(other, (tuple, list)) and len(other) >= 3:
            return other[0]*self[0] + other[1]*self[1] + other[2]*self[2]
        else:
            return NotImplemented

    def __sub__(self, other):
        """a - b: Vector subtraction"""
        if isinstance(other, (tuple, list)) and len(other) >= 3:
            return Vector(self[0]-other[0], self[1]-other[1], self[2]-other[2])
        else:
            return NotImplemented

    def __rsub__(self, other):
        """b - a: Vector subtraction"""
        if isinstance(other, (tuple, list)) and len(other) >= 3:
            return Vector(other[0]-self[0], other[1]-self[1], other[2]-self[2])
        else:
            return NotImplemented

    def __truediv__(self, other):
        """a / b: Scalar division"""
        if isinstance(other, (int, float)):
            return Vector(self[0]/other, self[1]/other, self[2]/other)
        else:
            return NotImplemented

    def __xor__(self, other):
        """a ^ b: Vector cross product"""
        if isinstance(other, (tuple, list)) and len(other) >= 3:
            return Vector(self[1]*other[2] - self[2]*other[1],
                          self[2]*other[0] - self[0]*other[2],
                          self[0]*other[1] - self[1]*other[0])
        else:
            return NotImplemented

    def __rxor__(self, other):
        """b ^ a: Vector cross product"""
        if isinstance(other, (tuple, list)) and len(other) >= 3:
            return Vector(other[1]*self[2] - other[2]*self[1],
                          other[2]*self[0] - other[0]*self[2],
                          other[0]*self[1] - other[1]*self[0])
        else:
            return NotImplemented

def find_fourth_vertex(row, distance1=2.45, distance2=1.53, distance3=2.50):
    # Use Vector type for the vertices
    p1 = Vector(row['x_coord_N'], row['y_coord_N'], row['z_coord_N'])
    p2 = Vector(row['x_coord_x'], row['y_coord_x'], row['z_coord_x']) # corresponds to the later step
    p3 = Vector(row['x_coord_C'], row['y_coord_C'], row['z_coord_C'])

    # Use float type for the distances
    r1 = float(distance1)
    r2 = float(distance2)
    r3 = float(distance3)

    u_axis = (p2 - p1).unit
    v_axis = (p3 - p1).perp(u_axis).unit
    w_axis = u_axis ^ v_axis

    u2 = (p2 - p1) | u_axis
    u3 = (p3 - p1) | u_axis
    v3 = (p3 - p1) | v_axis

    u = (r1*r1 - r2*r2 + u2*u2) / (2*u2)
    v = (r1*r1 - r3*r3 + u3*u3 + v3*v3 - 2*u*u3) / (2*v3)
    w = sqrt(r1*r1 - u*u - v*v)
    
    return p1 + u*u_axis + v*v_axis - w*w_axis
            # there is two correct mathematical solutions to the general problem but in our case only the second solution is correct - taking into account the constraints
            #(p1 + u*u_axis + v*v_axis + w*w_axis,
            #p1 + u*u_axis + v*v_axis - w*w_axis)
            
# calculate the Ca-Cb angle
def calculate_CaCb_angle(ligand_atoms,receptor_atoms, merged):
    
    # for CA → CB vector in receptor
    # convert to single letter 
    receptor_atoms_CACB = receptor_atoms[receptor_atoms.loc[:,"atom_name"].isin(["CA","CB"])].copy()
    receptor_atoms_CACB.loc[:,'receptor_amino_acid'] = receptor_atoms['residue_name'].map(amino_acids_dict)
    receptor_atoms_CACB.drop(["residue_name","chain_id","residue_number"], axis = 1, inplace = True)
    
    # for CA → CB vector in ligand
    ligand_atoms_CACB =  ligand_atoms[ligand_atoms.loc[:,"atom_name"].isin(["CA","CB"])].copy()
    ligand_atoms_CACB.drop(["chain_id"], axis = 1, inplace = True)
    ligand_CA = ligand_atoms_CACB[ligand_atoms_CACB["atom_name"] == "CA"]
    ligand_CB = ligand_atoms_CACB[ligand_atoms_CACB["atom_name"] == "CB"]

    ligand_atoms_CA_CB = pd.merge(ligand_CA, ligand_CB,
            on ="residue_number",
            how = "left")
    
    # calculate the vector in ligand
    # ligand_vector_angle → the CACB vector in ligand!! 
    ligand_atoms_CA_CB['ligand_vector_angle'] = ligand_atoms_CA_CB.apply(calculate_vector , axis = 1)
    ligand_residue_CA_CB = ligand_atoms_CA_CB[["residue_number", "ligand_vector_angle"]].copy()
         
    # Check if GLY in ligand
    if ligand_residue_CA_CB.isna().any(axis=1).any():
        print("Warning! GLY found in LIGAND sequence\n")    

    # calculate the vector in receptor
    receptor_CA = receptor_atoms_CACB[receptor_atoms_CACB["atom_name"] == "CA"]
    receptor_CB = receptor_atoms_CACB[receptor_atoms_CACB["atom_name"] == "CB"]

    receptor_atoms_CA_CB = pd.merge(receptor_CA, receptor_CB,
            on ="display_generic_number",
            how = "left")

    receptor_atoms_CA_CB['receptor_vector_angle'] = receptor_atoms_CA_CB.apply(calculate_vector , axis = 1)
    receptor_residue_CA_CB = receptor_atoms_CA_CB[["display_generic_number", "receptor_amino_acid_x", "receptor_vector_angle"]].copy()
    receptor_residue_CA_CB = receptor_residue_CA_CB.rename(columns={
        receptor_residue_CA_CB.columns[0]: "receptor_residue_generic_number",
        receptor_residue_CA_CB.columns[1]: "receptor_amino_acid"
    })
    receptor_residue_CA_CB.rename(columns={'receptor_amino_acid_x': 'receptor_amino_acid'}, inplace=True)

    # Check if GLY in RECEPTOR
    if receptor_residue_CA_CB.isna().any(axis=1).any():
        print("Warning! GLY found in RECEPTOR sequence\n")    
    
    # below deal with situations when Gly exist and calculate the Ca-Cb angle
    merged["residue_name"] = merged["residue_name"].map(amino_acids_dict)
    #  Filter rows in ligand/receptor
    ne = merged[merged["chain_id"].isin(['A', 'B', 'R'])]
    # Filter for GLY residues with 'N' and 'C' atoms and select relevant columns
    ne_N = ne[(ne["residue_name"] == "G") & (ne["atom_name"] == "N")].loc[:, ["atom_name", "chain_id","residue_number", "x_coord", "y_coord", "z_coord", "display_generic_number"]]
    ne_C = ne[(ne["residue_name"] == "G") & (ne["atom_name"] == "C")].loc[:, ["atom_name", "chain_id","residue_number", "x_coord", "y_coord", "z_coord", "display_generic_number"]]
    
    #  Merge the data frames for N and C atoms
    merged_ne_NC = pd.merge(ne_N, ne_C, on=["display_generic_number","chain_id", "residue_number"], how="left", suffixes=('_N', '_C'))

    CA_gly = merged[(merged["atom_name"] == "CA") & (merged["residue_name"] == "G")].copy()
    CA_gly.rename(columns={'x_coord': 'x_coord_x', 'y_coord': 'y_coord_x', 'z_coord': 'z_coord_x'}, inplace=True)
    
    # Calculate the coordinates of CB using the coordinates of C, CA, and N atoms
    merged_ne_final = pd.merge(CA_gly, merged_ne_NC, on=["display_generic_number","chain_id" ,"residue_number"], how="left")
    merged_ne_final["receptor_coordinates_CB"] = merged_ne_final.apply(find_fourth_vertex, axis=1)
    
    CB_gly = merged_ne_final.loc[:, ["display_generic_number", "chain_id", "residue_number", "residue_name", "receptor_coordinates_CB"]]
    CB_gly["atom_name"] = "CB"
    CB_gly[['x_coord_y', 'y_coord_y', 'z_coord_y']] = pd.DataFrame(CB_gly['receptor_coordinates_CB'].tolist(), index = CB_gly.index)
    CB_gly = CB_gly.drop('receptor_coordinates_CB', axis=1)
    
    all_CA = ne[(ne["chain_id"].isin(["A", "R", "B"])) & (ne["atom_name"] == "CA") & (ne["residue_name"] == "G") ]
    all_CA = all_CA[~((all_CA['chain_id'].isin(['A', 'R'])) & (all_CA['display_generic_number'].isnull()))]
    
    # integrate the coordinates of CA and CB atoms of GLY and calculate the CA-CB vector
    CA_CB_gly = pd.merge(all_CA, CB_gly, on=["chain_id", "residue_number" ,"display_generic_number"], how="left")
    CA_CB_gly.rename(columns={'x_coord': 'x_coord_x', 'y_coord': 'y_coord_x', 'z_coord': 'z_coord_x'}, inplace=True)
    CA_CB_gly['CA_CB_vector'] = CA_CB_gly.apply(calculate_vector, axis = 1)
    CA_CB_gly['amino_acid'] = "G"
    CA_CB_GLY = CA_CB_gly[["chain_id","display_generic_number", "residue_number","CA_CB_vector","amino_acid"]].copy()
    
    #  the CA-CB vector of the receptor's GLY residue   
    CA_CB_GLY_receptor = CA_CB_GLY[(CA_CB_GLY["chain_id"] == "A") | (CA_CB_GLY["chain_id"] == "R")].copy()
    CA_CB_GLY_receptor.rename(columns={'display_generic_number': 'receptor_residue_generic_number', 'amino_acid': 'receptor_amino_acid','CA_CB_vector':'receptor_vector_angle'}, inplace=True)
    CA_CB_GLY_receptor = CA_CB_GLY_receptor.drop(["chain_id", "residue_number"], axis=1)
    CaCb_vector_receptor = pd.concat([receptor_residue_CA_CB, CA_CB_GLY_receptor]).dropna().reset_index(drop=True)
    
    # the CA-CB vector of the ligand's GLY residue   
    CA_CB_GLY_ligand = CA_CB_GLY[(CA_CB_GLY["chain_id"] == "B")].copy()
    CA_CB_GLY_ligand.rename(columns={ 'amino_acid': 'receptor_amino_acid','CA_CB_vector':'ligand_vector_angle'}, inplace=True)
    CA_CB_GLY_ligand = CA_CB_GLY_ligand.drop(["chain_id", "display_generic_number","receptor_amino_acid"], axis=1)
    CaCb_vector_ligand = pd.concat([ligand_residue_CA_CB, CA_CB_GLY_ligand]).dropna().reset_index(drop=True)
    
    # CA_CB & CA_CB angle
    CaCb_vector_ligand['key'] = 1
    CaCb_vector_receptor['key'] = 1
    angle_CA_CB = pd.merge(CaCb_vector_receptor, CaCb_vector_ligand, on='key').drop('key', axis=1).copy()
    angle_CA_CB['Angle'] = angle_CA_CB.apply(calculate_angle, axis=1)
    angle_CA_CB.drop(["receptor_vector_angle","ligand_vector_angle"], axis = 1 ,inplace= True)
    
    return angle_CA_CB

# calculate the Toptwofrequency residue for each backbone
def prediction_ligands(Ca_Ca_distance, angle, reference_GPCR):
    # combine backbone result and angle result
    inference_backbone = pd.merge(Ca_Ca_distance, angle, 
                                  on=["receptor_residue_generic_number", "residue_number", "receptor_amino_acid"],
                                  how="left")
    
    # For each inference backbone residue, attempt to place each of the 20 natural amino acids in that position
    amino_acids = ['A', 'R', 'N', 'D', 'C', 'Q', 'E', 'G', 'H', 'I', 'L', 'K', 'M', 'F', 'P', 'S', 'T', 'W', 'Y', 'V']
    AA = pd.DataFrame(amino_acids, columns=['peptide_amino_acid'])
    inference_backbone["key"] = 1
    AA["key"] = 1
    peptide_residue_scores = pd.merge(inference_backbone, AA, on='key').drop('key', axis=1)

    final_df = pd.merge(peptide_residue_scores, reference_GPCR,
                        on=["receptor_residue_generic_number", "receptor_amino_acid", "peptide_amino_acid"],
                        how="left")

    # generate the score list table
    score_list = final_df[final_df['Distance_median'].notnull() & final_df['Angle_median'].notnull()] 
    
    # peptide_AA that fulfill the requirements for distance and angle
    peptide_AA_score = score_list.query('(Distance >= Distance_median - Distance_sd) & (Distance <= Distance_median + Distance_sd) & (Angle >= Angle_median - Angle_sd) & (Angle <= Angle_median + Angle_sd)')

    # add the frequencey and count together when RGN - residue_number - peptide_amnio_acid the same
    peptide_AA_score = peptide_AA_score.groupby(['residue_number', 'peptide_amino_acid']).agg(
        {
            'frequency(%)': 'sum',
            'count': 'sum'
        }
    ).reset_index()
    
    final_score = peptide_AA_score[["residue_number", "peptide_amino_acid", "frequency(%)", "count"]].sort_values(["residue_number", "frequency(%)"], ascending = [True,False])
    toptwofreqAA = final_score.groupby("residue_number").apply(lambda x: x.nlargest(2, "frequency(%)")).reset_index(drop=True)
    
    return toptwofreqAA

# generate predicted sequences based on Toptwofrequency
def generate_permutations(toptwofreqAA, backbone, group):
    starting_backbone_residue_number = min(group)
    ending_backbone_residue_number = max(group)
    
    # check for missing positions and add '.' for missing residues
    for i in range(starting_backbone_residue_number, ending_backbone_residue_number + 1):  
        if i not in toptwofreqAA['residue_number'].values:
            new_row = {'residue_number': i, 'peptide_amino_acid': '.', 'frequency(%)': 0, 'count': 0}
            toptwofreqAA = toptwofreqAA.append(new_row, ignore_index=True)

    # sort the DataFrame by residue_number and frequency
    toptwofreqAA = toptwofreqAA.sort_values(["residue_number", "frequency(%)"], ascending = [True,False]).copy()
    
    # filter residue number fulfill the group number
    toptwofreqAA = toptwofreqAA[(toptwofreqAA['residue_number'] >= starting_backbone_residue_number) & (toptwofreqAA['residue_number'] <= ending_backbone_residue_number)]
    
    # group by 'residue_number' to get peptide_amino_acid and corresponding frequencies/count
    grouped = toptwofreqAA.groupby('residue_number')
    position_data = {residue_number: group[['peptide_amino_acid', 'frequency(%)','count']].values.tolist() for residue_number, group in grouped}
    start_key = min(position_data.keys())

    # Step 1: generate unique combinations
    base_combination_chars = [chars[0][0] for chars in position_data.values()]
    unique_combinations = []

    for position, base_char in enumerate(base_combination_chars):
        current_combination = base_combination_chars.copy()
        if len(position_data[start_key + position]) > 1:
            current_combination[position] = position_data[start_key + position][1][0]
        combination_str = ''.join(current_combination)
        if combination_str not in unique_combinations:
            unique_combinations.append(combination_str)

    # Step 2: retrieve frequency for every position in unique combination
    def get_scores_for_combination(combination_str):
        char_scores = []
        for index, char in enumerate(combination_str):
            for char_info in position_data[start_key + index]:
                if char_info[0] == char:
                    char_scores.append(char_info[1])
                    break
        return char_scores

    scores_for_each_combination = {combination: get_scores_for_combination(combination) for combination in unique_combinations}

    # Step 3: calculate exact count sum for each combination
    integrated_results = []

    for combination, scores in scores_for_each_combination.items():
        exact_score_sum = sum(position_data[start_key + i][0][2] if position_data[start_key + i][0][0] == char else position_data[start_key + i][1][2] for i, char in enumerate(combination))
        integrated_results.append((combination, scores, exact_score_sum))

    result_df = pd.DataFrame(integrated_results, columns=['sequence', 'frequency_by_position','sum_count']) 

    # obtain the corresponding postion for the peptide backbone
    residue_num = list(range(starting_backbone_residue_number, ending_backbone_residue_number + 1))

    # convert the result to a dataFrame and rearrange columns
    result_df["residue_num"] = [residue_num] * len(result_df)
    result_df["struct_info"] = backbone

    # rearrangment
    col = result_df.pop('struct_info')
    result_df.insert(loc=0, column='struct_info', value=col)
    col_2 = result_df.pop('residue_num')
    result_df.insert(loc=1, column='residue_num', value=col_2)

    result_df["sum_frequency"] = result_df['frequency_by_position'].apply(sum)

    result_df = result_df.sort_values("sum_frequency", ascending=False).reset_index(drop=True)
    
    return result_df

# collect all possibile results for groups from one backbone
def merge_permutation_results(groups, toptwofreqAA, backbone):
    result_dfs = []  # list for storing the results of each group

    for group in groups:
        filtered_toptwofreqAA = toptwofreqAA.loc[toptwofreqAA['residue_number'].isin(group)].copy()      
        result_df = generate_permutations(filtered_toptwofreqAA, backbone, group)
        result_dfs.append(result_df)
    
    merged_df = pd.concat(result_dfs, axis=0).reset_index(drop=True)
    
    return merged_df

# find all groups from interaction information from one backbone
def find_interaction_part(s, min_length = 5):
    sequences = []  
    if not len(s): 
        print("No consecutive sequence found because the input array is empty.")
        return sequences  # Return an empty list

    current_sequence = [s[0]]  

    for i in range(1, len(s)):
        if s[i] == current_sequence[-1] + 1:
            current_sequence.append(s[i])  
        else:
            if len(current_sequence) >= min_length:
                sequences.append(current_sequence)  
            current_sequence = [s[i]] 
            
    if len(current_sequence) >= min_length:
        sequences.append(current_sequence)

    return sequences


# In[9]:


# 2.2 Implementation

from tqdm import tqdm
import os
import re
import requests
import pandas as pd
import numpy as np
from biopandas.pdb import PandasPdb
import requests
import pandas as pd

def list_files_in_directory(directory):
    entries = os.listdir(directory)
    files = [entry for entry in entries if os.path.isfile(os.path.join(directory, entry))]

    return files

directory_path = r"D:\GloriamGroup Dropbox\Sun Yifan\Yifan_peptide_project\Data\Experimental Structures Update"
file_list = list_files_in_directory(directory_path)

# original oprd_human_4RWD.pdb
backbones = [] # oprd_human_4RWD
receptors = [] # oprd_human
ligands = [] # ligand is from this PDB file: 4RWD
results_list = []
ligand_residues = {}

for item in file_list:
    parts = item.split('_')
    backbone = item[:-4]   
    receptor = parts[0] + '_' + parts[1]  
    ligand = parts[2].split('.')[0]  

    backbones.append(backbone)
    receptors.append(receptor)
    ligands.append(ligand)

for i in tqdm(range(len(file_list)), desc="Processing files"):
    backbone = backbones[i]
    entry_name = receptors[i]
    ligand_name = ligands[i]
    process_success = True  

    # attempt to read the PDB file and obtain the list of receptor residues
    try:
        inference_data = PandasPdb().read_pdb(f"D:\\GloriamGroup Dropbox\\Sun Yifan\\Yifan_peptide_project\\Data\\Experimental Structures Update\\{backbone}.pdb")
        # Check if there are any values in the 'insertion' column
        if (inference_data.df['ATOM']['insertion'].notnull() & (inference_data.df['ATOM']['insertion'] != '')).any():
            print('Warning: (Insertion Code) "insertion" column contains non-empty, non-null values!')
        # Check if there are any values in the 'insertion' column
        if (inference_data.df['HETATM']['insertion'].notnull() & (inference_data.df['HETATM']['insertion'] != '')).any():
            print('Warning : (Insertion Code) "insertion" column contains non-empty, non-null values!') 
       
        # only for nature residues
        df_1 = inference_data.df['ATOM']
        group_keys = ['atom_name', 'residue_name', 'chain_id', 'residue_number']
        sorted_data = df_1.sort_values(by=['occupancy', 'b_factor'], ascending=[False, True])
        df_1 = sorted_data.drop_duplicates(subset=group_keys, keep='first').copy()

        # also include unnature residues "X"
        df_2 = inference_data.df['HETATM']
        df_2.loc[:,"residue_name"] = 'X'
        sorted_data_2 = df_2.sort_values(by=['occupancy', 'b_factor'], ascending=[False, True])
        df_2 = sorted_data_2.drop_duplicates(subset=group_keys, keep='first').copy()

        inference_df = pd.concat([df_1, df_2])
        inference_df = inference_df[["atom_name", "residue_name", "chain_id", "residue_number", "x_coord", "y_coord", "z_coord"]].copy()

        
        # get the ligand sequence
        inference_df = inference_df.sort_values(by=["residue_number"])
        ligand = inference_df[(inference_df["chain_id"] == "B") & (inference_df["atom_name"] == "CA")]
        
        if ligand.empty:
            raise ValueError("No peptide column in the structural file")

        ligand = ligand.copy()
        ligand.loc[:, 'residue_name'] = ligand['residue_name'].map(amino_acids_dict)
        ligand_unique = ligand.drop_duplicates(subset=['residue_number'], keep='first')

        # 得到sequence和position
        ligand_sequence = ''.join(ligand_unique.loc[:,"residue_name"].tolist()) # 需要有顺序
        ligand_positions = np.array(ligand_unique["residue_number"].tolist())
        
        ligand_residues[backbone] = {
            'sequence': ligand_sequence,
            'positions': ligand_positions
        }
        
        # get the corresponding residue generic number
        url_2 = f"https://gpcrdb.org/services/residues/{entry_name}/"
        response_2 = requests.get(url_2)
        response_2.raise_for_status()  
        df_protein_residue = pd.DataFrame(response_2.json())
        df_protein_residue = df_protein_residue[["sequence_number", "display_generic_number"]]
        df_protein_residue = df_protein_residue.rename(columns={"sequence_number": "residue_number"})
    except Exception as e:
        print(f"Error occurred while initializing data for {backbone}, {entry_name}: {e}")
        process_success = False
        
    # attempt to obtain the interaction list and identify the appropriate groups
    if process_success:
        try:
            url = f"https://jimmy.gpcrdb.org/services/structure/{ligand_name}/peptideinteraction/"
            data_code = requests.get(url).json()
            peptide_residue_seq = pd.DataFrame(data_code)
            if 'peptide_residue_number' not in peptide_residue_seq.columns:
                raise ValueError("'No interaction data for the strcture in database")
            #  make predictions only for residue numbers that match the interaction requiremnt
            peptide_residue_seq_num = np.unique(peptide_residue_seq["peptide_residue_number"].tolist())
            groups = find_interaction_part(peptide_residue_seq_num)
            if not groups:
                raise ValueError("No consecutive sequences meeting the criteria were found")
        except Exception as e:
            print(f"Error occurred while finding groups for {backbone}, {ligand_name}: {e}")
            process_success = False

    # process the analysis and merge the results
    if process_success:
        try:
            ligand_atoms, receptor_atoms, receptor_atoms_CA, ligand_atoms_CA, merged = pre_treatment(inference_df, df_protein_residue)
            Ca_Ca_distance = calculate_CaCa_distance(receptor_atoms_CA, ligand_atoms_CA)
            angle = calculate_CaCb_angle(ligand_atoms, receptor_atoms, merged)
            toptwofreqAA = prediction_ligands(Ca_Ca_distance, angle, reference_GPCR)
            result = merge_permutation_results(groups, toptwofreqAA, backbone)
            results_list.append(result)

        except Exception as e:
            print(f"An error occurred during the analysis for {backbone}, {entry_name}: {e}")
         
result_df = pd.concat(results_list, ignore_index=True)


# In[13]:


# 2.3 Scoring and evaluation

# obtain the actual amino acid sequence at the corresponding positions

def extract_truth_seq(ligand_residues, query_key, query_positions):
    if query_key not in ligand_residues:
        return "Query key not found in ligand residues."
    
    query_amino_acids = ""
    for pos in query_positions:
        try:
            index = np.where(np.array(ligand_residues[query_key]["positions"]) == pos)[0][0]
            query_amino_acids += ligand_residues[query_key]["sequence"][index]
        except IndexError:
#             print(f"Warning: Position {pos} not found in {query_key}. Adding '-' as placeholder.")
            query_amino_acids += "-"
        except Exception as e:
            return f"An error occurred: {str(e)}"
    
    return query_amino_acids

result_df['origin_seq'] = result_df.apply(lambda row: extract_truth_seq(ligand_residues, row['struct_info'], row['residue_num']), axis=1)

# scoring method
# calculating the median and standard deviation
result_df['median_by_position'] = result_df['frequency_by_position'].apply(np.median)
result_df['std_by_position'] = result_df['frequency_by_position'].apply(np.std)
result_df['cv_by_position'] = result_df['std_by_position'] / result_df['median_by_position']

# length of valid characters in the sequence 
def filter_condition(row):
    valid_sequence = row['sequence'].replace('.', '')  # remove "."
    valid_origin_seq = row['origin_seq'].replace('-', '')  # remove "-"
    return len(valid_sequence) >= 5 and len(valid_origin_seq) >= 5

filtered_df = result_df[result_df.apply(filter_condition, axis=1)]
result_df = filtered_df.copy()

# evaluation
def compare_sequences(seq1, seq2, substitution_matrix):
    try:
        if len(seq1) != len(seq2):
            raise ValueError("Sequence length NOT EQUAL")

        clusters = {
            "STAGP": set("STAGP"),
            "DEQN": set("DEQN"),
            "HRK": set("HRK"),
            "MILV": set("MILV"),
            "WFY": set("WFY"),
            "LFI": set("LFI"),
            "FHY": set("FHY")
        }

        identity_score = 0
        similarity_score = 0
        Blosum62_score = 0
        
        for a, b in zip(seq1, seq2):
            if a in ['-', '.'] or b in ['-', '.','X']:    # b indicates seq2 
                continue  

            if a == b:
                identity_score += 1
                similarity_score += 1
            else:
                for cluster in clusters.values():
                    if a in cluster and b in cluster:
                        similarity_score += 1
                        break

            if (a, b) in substitution_matrix:
                Blosum62_score += substitution_matrix[a, b]
            else:
                raise KeyError(f"Substitution matrix missing value for pair: ({a}, {b})")

        seq_length = len(seq2.replace("X", ""))  # Using the length WITHOUT X of seq2 for calculation

        return (
            identity_score / seq_length,
            similarity_score / seq_length,
            Blosum62_score
        )


    except ValueError as ve:
        print(ve)
        return None, None, None
    except KeyError as ke:
        print(f"Error: No key in substitution matrix {ke}")
        return None, None, None

# modify the Blosum62 matrix   
from Bio.Align import substitution_matrices
substitution_matrix = substitution_matrices.load('BLOSUM62') 
# change the diagonal value to 5 for specified amino acids
amino_acids = 'CSTPAGNDEQHRKMILVFYW'
for aa in amino_acids:
    substitution_matrix[aa, aa] = 5
    
# scores for 'X' against all amino acids to 0
for aa in amino_acids:
    substitution_matrix['X', aa] = 0
    substitution_matrix[aa, 'X'] = 0

final_tab = result_df.copy()
final_tab[['identity', 'similarity',"Blosum62_m"]] = final_tab.apply(
    lambda x: compare_sequences(x['sequence'], x['origin_seq'],substitution_matrix), 
    axis=1, result_type="expand"
)

final_tab["freq_perAA"] = final_tab.apply(
    lambda x: x["sum_frequency"] / len(x["sequence"]) if len(x["sequence"]) > 0 else 0, axis=1)
final_tab["count_perAA"] = final_tab.apply(
    lambda x: x["sum_count"] / len(x["sequence"]) if len(x["sequence"]) > 0 else 0, axis=1)

validation_table = final_tab[~final_tab['origin_seq'].str.contains('X')].copy()
validation_table


# # Part3: Prediction 

# In[15]:


# 3.1 Require functions

def generate_permutations(toptwofreqAA, backbone):
    # determine the starting and ending residue numbers from the input data
    starting_backbone_residue_number = toptwofreqAA['residue_number'].min()
    ending_backbone_residue_number = toptwofreqAA['residue_number'].max()
    
    # check for missing positions and add '.' for missing residues
    for i in range(starting_backbone_residue_number, ending_backbone_residue_number + 1):  
        if i not in toptwofreqAA['residue_number'].values:
            new_row = {'residue_number': i, 'peptide_amino_acid': '.', 'frequency(%)': 0, 'count': 0}
            toptwofreqAA = toptwofreqAA.append(new_row, ignore_index=True)

    toptwofreqAA = toptwofreqAA.sort_values(["residue_number", "frequency(%)"], ascending=[True, False]).copy()

    grouped = toptwofreqAA.groupby('residue_number')
    position_data = {residue_number: group[['peptide_amino_acid', 'frequency(%)', 'count']].values.tolist() for residue_number, group in grouped}

    start_key = min(position_data.keys())

    # Step 1: generate unique combinations
    base_combination_chars = [chars[0][0] for chars in position_data.values()]
    unique_combinations = []

    for position, base_char in enumerate(base_combination_chars):
        current_combination = base_combination_chars.copy()
        if len(position_data[start_key + position]) > 1:
            current_combination[position] = position_data[start_key + position][1][0]
        combination_str = ''.join(current_combination)
        if combination_str not in unique_combinations:
            unique_combinations.append(combination_str)

    # Step 2: retrieve frequency for every position in unique combination
    def get_scores_for_combination(combination_str):
        char_scores = []
        for index, char in enumerate(combination_str):
            for char_info in position_data[start_key + index]:
                if char_info[0] == char:
                    char_scores.append(char_info[1])
                    break
        return char_scores

    scores_for_each_combination = {combination: get_scores_for_combination(combination) for combination in unique_combinations}

    # Step 3: calculate exact count sum for each combination
    integrated_results = []

    for combination, scores in scores_for_each_combination.items():
        exact_score_sum = sum(position_data[start_key + i][0][2] if position_data[start_key + i][0][0] == char else position_data[start_key + i][1][2] for i, char in enumerate(combination))
        integrated_results.append((combination, scores, exact_score_sum))

    # creating a DataFrame from the combinations with frequency and count sum
    result_df = pd.DataFrame(integrated_results, columns=['sequence', 'frequency_by_position', 'sum_count'])
    residue_num = list(range(starting_backbone_residue_number, ending_backbone_residue_number + 1))

    # convert the result to a DataFrame and rearrange columns
    result_df["residue_num"] = [residue_num] * len(result_df)
    result_df["struct_info"] = backbone

    col = result_df.pop('struct_info')
    result_df.insert(loc=0, column='struct_info', value=col)
    col_2 = result_df.pop('residue_num')
    result_df.insert(loc=1, column='residue_num', value=col_2)
    result_df["sum_frequency"] = result_df['frequency_by_position'].apply(sum)
    result_df = result_df.sort_values("sum_frequency", ascending=False).reset_index(drop=True)
    
    return result_df

def merge_permutation_results(toptwofreqAA, backbone):
    merged_df = generate_permutations(toptwofreqAA, backbone)
    
    return merged_df


# In[16]:


# 3.2 Implementation
file_list = ['7F1T_aa2ar_human_1.pdb',
 '7WVY_gpr62_human_1.pdb',
 '7F1T_aa2br_human_1.pdb',
 '7WVY_gpr52_human_1.pdb',
 '7XBX_gpr61_human_1.pdb',
 '7XBX_5ht7r_human_1.pdb',
 '7XBX_gp161_human_1.pdb',
 '7W0O_gasr_human_2.pdb',
 '8F7X_gpr78_human_1.pdb',
 '7F8V_gp107_human_1.pdb',
 '7XBX_gpr26_human_1.pdb',
 '7W0O_5ht2a_human_1.pdb',
 '8DWC_par2_human_1.pdb',
 '7VLA_gp135_human_2.pdb',
 '7VLA_gp135_human_2.pdb',
 '7P02_mas1l_human_1.pdb',
 '6JOD_aa1r_human_1.pdb',
 '7F1T_aa1r_human_1.pdb',
 '7F53_c5ar2_human_1.pdb',
 '8F7Q_o51d1_human_1.pdb']


# original oprd_human_4RWD.pdb
backbones = [] # oprd_human_4RWD
receptors = [] # oprd_human

for item in file_list:
    parts = item.split('_')
    backbone = item[:-4]   
    receptor = item[5:-6]   
    
    backbones.append(backbone)
    receptors.append(receptor)


from biopandas.pdb import PandasPdb
import requests
import pandas as pd

results_list = []
ligand_residues = {}

for i in tqdm(range(len(file_list)), desc="Processing files"):
    backbone = backbones[i]
    entry_name = receptors[i]
    process_success = True  


    try:

        inference_data = PandasPdb().read_pdb(f"D:\\GloriamGroup Dropbox\\Sun Yifan\\Yifan_peptide_project\\Data\\Superimposed_structures_6Dec_ALA\\{backbone}.pdb")
        # check if there are any values in the 'insertion' column
        if (inference_data.df['ATOM']['insertion'].notnull() & (inference_data.df['ATOM']['insertion'] != '')).any():
            print('Warning: (Insertion Code) "insertion" column contains non-empty, non-null values!')

        # check if there are any values in the 'insertion' column
        if (inference_data.df['HETATM']['insertion'].notnull() & (inference_data.df['HETATM']['insertion'] != '')).any():
            print('Warning : (Insertion Code) "insertion" column contains non-empty, non-null values!') 
        # only for nature residues
        df_1 = inference_data.df['ATOM']
        group_keys = ['atom_name', 'residue_name', 'chain_id', 'residue_number']
        sorted_data = df_1.sort_values(by=['occupancy', 'b_factor'], ascending=[False, True])
        df_1 = sorted_data.drop_duplicates(subset=group_keys, keep='first').copy()

        # also include unnature residues "X"
        df_2 = inference_data.df['HETATM']
        df_2.loc[:,"residue_name"] = 'X'
        sorted_data_2 = df_2.sort_values(by=['occupancy', 'b_factor'], ascending=[False, True])
        df_2 = sorted_data_2.drop_duplicates(subset=group_keys, keep='first').copy()

        inference_df = pd.concat([df_1, df_2])
        inference_df = inference_df[["atom_name", "residue_name", "chain_id", "residue_number", "x_coord", "y_coord", "z_coord"]].copy()
    
        inference_df = inference_df.sort_values(by=["residue_number"])
        ligand = inference_df[(inference_df["chain_id"] == "B") & (inference_df["atom_name"] == "CA")]
        if ligand.empty:
            raise ValueError("No peptide column in the structural file")
        
        ligand = ligand.copy()
        ligand.loc[:, 'residue_name'] = ligand['residue_name'].map(amino_acids_dict)
        ligand_unique = ligand.drop_duplicates(subset=['residue_number'], keep='first')

        ligand_sequence = ''.join(ligand_unique.loc[:,"residue_name"].tolist()) 
        ligand_positions = np.array(ligand_unique["residue_number"].tolist())
        
        ligand_residues[backbone] = {
            'sequence': ligand_sequence,
            'positions': ligand_positions
        }
        
        url_2 = f"https://gpcrdb.org/services/residues/{entry_name}/"
        response_2 = requests.get(url_2)
        response_2.raise_for_status()  
        df_protein_residue = pd.DataFrame(response_2.json())
        df_protein_residue = df_protein_residue[["sequence_number", "display_generic_number"]]
        df_protein_residue = df_protein_residue.rename(columns={"sequence_number": "residue_number"})
    except Exception as e:
        print(f"Error occurred while initializing data for {backbone}, {entry_name}: {e}")
        process_success = False
        
    if process_success:
        try:
            ligand_atoms, receptor_atoms, receptor_atoms_CA, ligand_atoms_CA, merged = pre_treatment(inference_df, df_protein_residue)
            Ca_Ca_distance = calculate_CaCa_distance(receptor_atoms_CA, ligand_atoms_CA)
            angle = calculate_CaCb_angle(ligand_atoms, receptor_atoms, merged)
            toptwofreqAA = prediction_ligands(Ca_Ca_distance, angle, reference_GPCR)
            if toptwofreqAA.empty:
                    raise ValueError("All predicted residue out of the range requirement for CA and CB")
            
            result = merge_permutation_results(toptwofreqAA, backbone)
            results_list.append(result)

        except Exception as e:
            print(f"An error occurred during the analysis for {backbone}, {entry_name}: {e}")

result_df = pd.concat(results_list, ignore_index=True)


# In[20]:


# 3.3 Scring
def extract_truth_seq(ligand_residues, query_key, query_positions):
    if query_key not in ligand_residues:
        return "Query key not found in ligand residues."
    
    query_amino_acids = ""
    for pos in query_positions:
        try:
            index = np.where(np.array(ligand_residues[query_key]["positions"]) == pos)[0][0]
            query_amino_acids += ligand_residues[query_key]["sequence"][index]
        except IndexError:
            query_amino_acids += "-"
        except Exception as e:
            return f"An error occurred: {str(e)}"
    
    return query_amino_acids

result_df['origin_seq'] = result_df.apply(lambda row: extract_truth_seq(ligand_residues, row['struct_info'], row['residue_num']), axis=1)

result_df['median_by_position'] = result_df['frequency_by_position'].apply(np.median)
result_df['std_by_position'] = result_df['frequency_by_position'].apply(np.std)
result_df['cv_by_position'] = result_df['std_by_position'] / result_df['median_by_position']

def filter_condition(row):
    valid_sequence = row['sequence'].replace('.', '')  
    valid_origin_seq = row['origin_seq'].replace('-', '')  
    return len(valid_sequence) >= 5 and len(valid_origin_seq) >= 5

# filter
filtered_df = result_df[result_df.apply(filter_condition, axis=1)]
result_df = filtered_df.copy()

result_df["freq_perAA"] = result_df.apply(
    lambda x: x["sum_frequency"] / len(x["sequence"]) if len(x["sequence"]) > 0 else 0, axis=1)

result_df["count_perAA"] = result_df.apply(
    lambda x: x["sum_count"] / len(x["sequence"]) if len(x["sequence"]) > 0 else 0, axis=1)

prediction_table = result_df.copy()
prediction_table


# In[90]:





# In[ ]:




