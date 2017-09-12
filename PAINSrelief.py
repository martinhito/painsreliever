from __future__ import print_function
import rdkit
from rdkit import Chem
from rdkit.Chem import rdmolops
from rdkit import DataStructs
from rdkit.Chem.Fingerprints import FingerprintMols
from rdkit.Chem import MACCSkeys
import csv



########## Populating arrays from CSV files ########## 

painsData = []	#tuples of data about everyPAINS alert taken from master CSV

with open('C:/Users/marti/Downloads/PAINS_master.csv') as painsAlertCSV:
	reader = csv.DictReader(painsAlertCSV, delimiter=',')
	for row in reader:
		tempDataTuple =(row['SMARTS'], row['\xef\xbb\xbfALERT NAME'], row['NPubChem'], row['NDCM'], row['Luciferase'], row['\xce\xb2-lactamase'], row['Fluorescence'], row['All Assays']) #SMARTS string for each alert is at index 0
		painsData.append(tempDataTuple)
		
compoundData = [] #array of tuples, everything about compounds, to be used for MACCSkeys analysis

with open('C:/Users/marti/Downloads/PAINS_compounds_noDuplicates.csv') as compoundsCSV:
	reader = csv.DictReader(compoundsCSV, delimiter=',')
	for row in reader:
		tempDataTuple =((row['\xef\xbb\xbfCID'], row['SMILES'], row['AllAssays_Active'], row['AllAssays_Total'], row['LuciferaseAssays_Active'], row['LuciferaseAssays_Total'], row['BetaLactamaseAssays_Active'], row['BetaLactamaseAssays_Total'], row['FluorescenceAssays_Active'], row['FluorescenceAssays_Total'], row['PAINS_A'], row['PAINS_B'], row['PAINS_C'])) #3 PAINS alerts b/c a compound can have more than 1
		compoundData.append(tempDataTuple)


		
########## Flagging PAINS Alerts on query compound ##########

matchIndices = []	#indices of flagged alert on molecule (for highlighting)
flaggedAlerts = []	#information about each FLAGGED PAINS alert, taken from painsData[]

querySmiles = 'Fc1c(c(ccc1)Cl)C=C2SC(=S)N(C2=O)C(C)c3ccccc3' #user-inputted SMILES string
queryMol = Chem.MolFromSmiles(querySmiles) #converted to mol

for alert in painsData:
	painsAlert = Chem.MolFromSmarts(alert[0]) #converting alert SMARTS string to mol (SMARTS located at index 0 of alert tuple)
	
	if Chem.AddHs(queryMol).HasSubstructMatch(painsAlert): #search with explicit hydrogens added 
		matchIndices.append(Chem.AddHs(queryMol).GetSubstructMatch(painsAlert)) #adds indices of PAINS alert substructure on query molecule 
		flaggedAlerts.append(alert)
		
	elif queryMol.HasSubstructMatch(painsAlert): #search without explicit hydrogens added
		flaggedAlerts.append(alert)
		matchIndices.append(queryMol.GetSubstructMatch(painsAlert))
		
for match in flaggedAlerts:
	print(match)
	
	

########## Getting Compounds with same alert ##########

sameAlertCompounds = [] #all compounds with same alert as query molecule

for painsAlert in flaggedAlerts:
	for row in compoundData: #10,11,12 are where the PAINS alerts are in the tuples of flaggedAlerts
		if row[10] == painsAlert[1]: #painsAlert[1] is the alert name
			sameAlertCompounds.append(row)
		if row[11] == painsAlert[1]: 
			sameAlertCompounds.append(row)
		if row[12] == painsAlert[1]: 
			sameAlertCompounds.append(row)
	
	

	
########## MACCS Fingerprint analysis ########## 

highestTc = 0		#Tc is a measure of how similar the compound is to the query molecule
mostSimilarMol = 0	#compound with highest Tc value
similarMols = []	#all compounds with Tc > 0.9, not including the mostSimilarMol

queryFPS = MACCSkeys.GenMACCSKeys(Chem.MolFromSmiles(querySmiles)) #FPS=fingerprints which is how MACCSkeys analysis works

for compound in sameAlertCompounds:
	compoundFPS = MACCSkeys.GenMACCSKeys(Chem.MolFromSmiles(compound[1])) #1 is index of the SMILES string in the tuple 'compound' in sameAlertCompounds
	if DataStructs.FingerprintSimilarity(compoundFPS, queryFPS) > highestTc: #finding mostSimilarMol
		highestTc = DataStructs.FingerprintSimilarity(compoundFPS, queryFPS)
		mostSimilarMol = compound
	elif DataStructs.FingerprintSimilarity(compoundFPS, queryFPS) > 0.9: #Tc value > 0.9
		similarMols.append(compound)
similarMols.append(mostSimilarMol) #puts most similar mol at end of list

print(mostSimilarMol)
for compound in similarMols:
	print(compound)

	


					