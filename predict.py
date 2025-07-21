from esm.models.esmc import ESMC
import pandas as pd
from glycowork.ml.inference import get_lectin_preds, get_esmc_representations
from glycowork.ml.models import prep_model
import pandas as pd
from Bio import SeqIO
import argparse
import yaml

class BindingPredictor():
    def __init__(self, config):
        glycan_path = config["glycan_path"]
        protein_path = config["protein_path"]
        fasta_path = config["fasta_path"]

        self.glycan_df = pd.read_csv(glycan_path, sep="\t")
        self.protein_df = pd.read_csv(protein_path, sep="\t", index_col=0)
        self.protein_fasta = self.load_protein_fasta(fasta_path)
        
        self.esm_model = ESMC.from_pretrained("esmc_300m")
        self.leor = prep_model("LectinOracle", 1, trained = True)
        
        
    def process_iupacs(self, glycans):
        for count, iupac in enumerate(glycans):     
            #remove spacer
            if iupac.count('(') > iupac.count(')'):
                spacer_pos = iupac.rfind('(')
                iupac = iupac[:spacer_pos]
                    
            # format sulfate/phosphate groups
            while 'O' in iupac:
                oxygen_idx = iupac.find('O')
                ox_modification_idx = oxygen_idx + 1
                ox_modification = iupac[ox_modification_idx]
                
                end_idx = ox_modification_idx + 1
                while(end_idx < len(iupac) and iupac[end_idx].isnumeric()):
                    end_idx += 1
                num_groups = int(iupac[ox_modification_idx + 1: end_idx]) if ox_modification_idx + 1 < end_idx else 1

                positions = []
                start_idx = oxygen_idx
                for _ in range(num_groups):
                    start_idx -= 1
                    position_end = start_idx
                    while iupac[start_idx].isnumeric():
                        start_idx -= 1
                    positions.insert(0, iupac[start_idx + 1: position_end + 1])
                start_idx += 1
                replacement = "".join(f"{position}{ox_modification}" for position in positions)
                iupac = iupac.replace(iupac[start_idx : end_idx], replacement)
                
            # remove MDPLys spacer
            iupac = iupac.replace('MDPLys', '')
                
            glycans[count] = iupac
        return glycans

    def load_protein_fasta(self, fasta_file):
        # list to collect parsed records
        records_list = []

        for record in SeqIO.parse(fasta_file, "fasta"):
            header_parts = record.description.split(" ")

            uniprot_id = header_parts[0].split("|")[1]
            entry_name = header_parts[0].split("|")[2]

            # Get positions of OS= (organism species start)
            try:
                os_index = next(i for i, s in enumerate(header_parts) if s.startswith("OS="))
            except StopIteration:
                os_index = None

            # Protein name is from after entry_name (index 0) up to OS=
            protein_name = " ".join(header_parts[1:os_index])

            # Extract fields safely if present
            os = next((part.split("=")[1] for part in header_parts if part.startswith("OS=")), None)
            ox = next((part.split("=")[1] for part in header_parts if part.startswith("OX=")), None)
            gn = next((part.split("=")[1] for part in header_parts if part.startswith("GN=")), None)
            pe = next((part.split("=")[1] for part in header_parts if part.startswith("PE=")), None)
            sv = next((part.split("=")[1] for part in header_parts if part.startswith("SV=")), None)

            sequence = str(record.seq)

            records_list.append({
                "Accession": uniprot_id,
                "EntryName": entry_name,
                "ProteinName": protein_name,
                "Organism": os,
                "TaxonomyID": ox,
                "Gene": gn,
                "ProteinExistence": pe,
                "SequenceVersion": sv,
                "Sequence": sequence
            })

        return pd.DataFrame(records_list)
    
    def get_sequences(self):
        self.protein_df["Sequence"] = pd.NA 
        for i, row in self.protein_df.iterrows():
            protein_name = row["Name"]
            uniprot_id = row["UniProt ID"]
            if pd.isna(uniprot_id):
                continue
            
            if "GenBank" in uniprot_id:
                uniprot_id = uniprot_id.replace('GenBank', '').strip()
            if "/" in uniprot_id:
                uniprot_id = uniprot_id.split('/')[0].strip()
            if "(" in uniprot_id:
                uniprot_id = uniprot_id.split('(')[0].strip()
                
            matching_rows = self.protein_fasta.loc[self.protein_fasta['Accession'] == uniprot_id].reset_index()
            # if uniprot_id doesnt match Acession value, try EntryName value
            if len(matching_rows) != 1:
                matching_rows = self.protein_fasta.loc[self.protein_fasta['EntryName'] == uniprot_id].reset_index()
                
            # print(matching_rows)
            # print(self.protein_fasta)
            if len(matching_rows) != 1:
                print(f"Could not find sequence for protein '{protein_name}' with UniProtID '{uniprot_id}'")
            elif len(matching_rows) == 1:
                sequence = matching_rows["Sequence"][0]
                self.protein_df.at[i, 'Sequence'] = sequence
                
        self.protein_df = self.protein_df.dropna(subset=["Sequence"])

    def get_pred_df(self, protein_df, glycan_df):
        leor_preds = {"glycan": [], "protein": [], "rfu": []}
        for i, prot_row in protein_df.iterrows():
            sequence = prot_row["Sequence"]
            protein = prot_row["Acronym"]
            rep = get_esmc_representations([sequence], self.esm_model)
            preds = get_lectin_preds(sequence, glycan_df["Clean IUPAC"], self.leor, rep, sort=False)
            for j, pred_row in preds.iterrows():
                glycan = glycan_df["Name"][j]
                rfu = pred_row["pred"]

                leor_preds["glycan"].append(glycan)
                leor_preds["protein"].append(protein)
                leor_preds["rfu"].append(rfu)
                
        return pd.DataFrame(leor_preds)
    
    def predict(self, save_path):
        self.glycan_df["Clean IUPAC"] = self.process_iupacs(self.glycan_df["IUPAC"])
        
        self.protein_df = self.protein_df[["Protein Name", "Full protein name", "UniProt ID, etc."]]
        self.protein_df = self.protein_df.rename(columns={
            "UniProt ID, etc.": "UniProt ID", 
            "Full protein name": "Name", 
            "Protein Name": "Acronym"
        })
        self.get_sequences()
        
        print(f"\nPredicting the interactions of {len(self.glycan_df)} glycans with {len(self.protein_df)} lectins...")
        leor_preds_df = self.get_pred_df(self.protein_df, self.glycan_df)
        leor_preds_df.to_csv(save_path, index=0)
        print(f"\nInteraction predictions saved to {save_path}")

def main():
    parser = argparse.ArgumentParser(description="Predict protein-glycan interactions")
    parser.add_argument("--config", type=str, required=True, help="Path to the YAML configuration file")
    cli_args = parser.parse_args()
    
    with open(cli_args.config, 'r') as file:
        config = yaml.safe_load(file)
    
    predictor = BindingPredictor(config)
    predictor.predict(config["save_path"])
        
if __name__ == "__main__":
    main()


