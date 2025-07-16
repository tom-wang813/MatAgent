from rdkit import Chem, DataStructs
from rdkit.Chem import Descriptors, Crippen, Lipinski, AllChem
import json
from pydantic import BaseModel, Field
from typing import List, Literal

from mcp_server.core.registry import register_basic_tool
from mcp_server.core.schemas import BasicToolMetadata

def _get_mol(smiles_string: str):
    """Helper function to convert SMILES to a Mol object, with error handling."""
    mol = Chem.MolFromSmiles(smiles_string)
    if mol is None:
        raise ValueError(f"Invalid SMILES string provided: {smiles_string}")
    return mol

class CalculateAllDescriptorsInput(BaseModel):
    smiles_string: str = Field(..., description="The SMILES representation of the molecule.")

# --- Tool: Calculate All Descriptors ---
@register_basic_tool(
    metadata=BasicToolMetadata(
        name="calculate_all_descriptors",
        description="Calculates a comprehensive set of RDKit molecular descriptors for a given SMILES string.",
        tags=["rdkit", "chemistry", "descriptor", "comprehensive"],
        author="RDKit Wrapper",
        parameters=CalculateAllDescriptorsInput.model_json_schema()
    )
)
def calculate_all_descriptors(input: CalculateAllDescriptorsInput) -> dict:
    mol = _get_mol(input.smiles_string)
    descriptors = {
        "MolWt": Descriptors.MolWt(mol),
        "LogP": Crippen.MolLogP(mol),
        "TPSA": Descriptors.TPSA(mol),
        "NumRotatableBonds": Descriptors.NumRotatableBonds(mol),
        "NumHDonors": Lipinski.NumHDonors(mol),
        "NumHAcceptors": Lipinski.NumHAcceptors(mol),
        "NumHeavyAtoms": Descriptors.HeavyAtomCount(mol),
        "NumAromaticRings": Descriptors.NumAromaticRings(mol),
        "NumAliphaticRings": Descriptors.NumAliphaticRings(mol),
        "NumSaturatedRings": Descriptors.NumSaturatedRings(mol),
        "NumHeteroatoms": Descriptors.NumHeteroatoms(mol),
        "FractionCSP3": Descriptors.FractionCSP3(mol),
        "RingCount": Descriptors.RingCount(mol),
        "ExactMolWt": Descriptors.ExactMolWt(mol),
        "NumValenceElectrons": Descriptors.NumValenceElectrons(mol),
        "NumRadicalElectrons": Descriptors.NumRadicalElectrons(mol),
    }
    return descriptors

class GenerateDescriptorReportInput(BaseModel):
    smiles_string: str = Field(..., description="The SMILES representation of the molecule.")
    output_format: Literal["csv"] = Field("csv", description="The desired output format (currently only csv is supported).")

# --- Tool: Generate Descriptor Report ---
@register_basic_tool(
    metadata=BasicToolMetadata(
        name="generate_descriptor_report",
        description="Generates a CSV report of molecular descriptors for a given SMILES string.",
        tags=["rdkit", "chemistry", "descriptor", "report", "file-generation"],
        author="RDKit Wrapper",
        parameters=GenerateDescriptorReportInput.model_json_schema()
    )
)
def generate_descriptor_report(input: GenerateDescriptorReportInput) -> str:
    if input.output_format != "csv":
        raise ValueError("Unsupported output format. Only 'csv' is supported.")

    descriptors = calculate_all_descriptors(CalculateAllDescriptorsInput(smiles_string=input.smiles_string))

    # Prepare CSV content
    header = ",".join(descriptors.keys())
    values = ",".join(map(str, descriptors.values()))
    csv_content = f"{header}\n{values}"

    return csv_content

class CompareMoleculesReportInput(BaseModel):
    smiles_list: List[str] = Field(..., description="A list of SMILES representations of the molecules to compare.")
    properties: List[str] = Field(..., description="A list of descriptor names to include in the report (e.g., [\"MolWt\", \"LogP\"]).")
    output_format: Literal["csv"] = Field("csv", description="The desired output format (currently only csv is supported).")

# --- Tool: Compare Molecules Report ---
@register_basic_tool(
    metadata=BasicToolMetadata(
        name="compare_molecules_report",
        description="Compares a list of molecules based on specified RDKit descriptors and generates a CSV report.",
        tags=["rdkit", "chemistry", "descriptor", "comparison", "report", "file-generation"],
        author="RDKit Wrapper",
        parameters=CompareMoleculesReportInput.model_json_schema()
    )
)
def compare_molecules_report(input: CompareMoleculesReportInput) -> str:
    if input.output_format != "csv":
        raise ValueError("Unsupported output format. Only 'csv' is supported.")

    all_data = []
    for smiles_string in input.smiles_list:
        try:
            descriptors = calculate_all_descriptors(CalculateAllDescriptorsInput(smiles_string=smiles_string))
            row_data = {"SMILES": smiles_string}
            for prop in input.properties:
                row_data[prop] = descriptors.get(prop, "N/A") # Get property or N/A if not found
            all_data.append(row_data)
        except ValueError:
            # Handle invalid SMILES by adding a row with N/A for properties
            row_data = {"SMILES": smiles_string}
            for prop in input.properties:
                row_data[prop] = "Invalid SMILES"
            all_data.append(row_data)

    if not all_data:
        return ""

    # Prepare CSV content
    header = ["SMILES"] + properties
    csv_lines = [",".join(header)]
    for row in all_data:
        values = [str(row.get(h, "N/A")) for h in header]
        csv_lines.append(",".join(values))

    return "\n".join(csv_lines)

class CalculateMolecularWeightInput(BaseModel):
    smiles_string: str = Field(..., description="The SMILES representation of the molecule.")

# --- Tool: Calculate Molecular Weight ---
@register_basic_tool(
    metadata=BasicToolMetadata(
        name="calculate_molecular_weight",
        description="Calculates the molecular weight (g/mol) of a molecule given its SMILES string.",
        tags=["rdkit", "chemistry", "descriptor"],
        author="RDKit Wrapper",
        parameters=CalculateMolecularWeightInput.model_json_schema()
    )
)
def calculate_molecular_weight(input: CalculateMolecularWeightInput) -> float:
    mol = _get_mol(input.smiles_string)
    return Descriptors.MolWt(mol)

class CalculateLogPInput(BaseModel):
    smiles_string: str = Field(..., description="The SMILES representation of the molecule.")

# --- Tool: Calculate LogP ---
@register_basic_tool(
    metadata=BasicToolMetadata(
        name="calculate_logp",
        description="Calculates the octanol-water partition coefficient (LogP) of a molecule, a measure of its lipophilicity.",
        tags=["rdkit", "chemistry", "descriptor", "lipophilicity"],
        author="RDKit Wrapper",
        parameters=CalculateLogPInput.model_json_schema()
    )
)
def calculate_logp(input: CalculateLogPInput) -> float:
    mol = _get_mol(input.smiles_string)
    return Crippen.MolLogP(mol)

class CalculateTPSAInput(BaseModel):
    smiles_string: str = Field(..., description="The SMILES representation of the molecule.")

# --- Tool: Calculate TPSA ---
@register_basic_tool(
    metadata=BasicToolMetadata(
        name="calculate_tpsa",
        description="Calculates the Topological Polar Surface Area (TPSA) of a molecule.",
        tags=["rdkit", "chemistry", "descriptor", "polarity"],
        author="RDKit Wrapper",
        parameters=CalculateTPSAInput.model_json_schema()
    )
)
def calculate_tpsa(input: CalculateTPSAInput) -> float:
    mol = _get_mol(input.smiles_string)
    return Descriptors.TPSA(mol)

class CountRotatableBondsInput(BaseModel):
    smiles_string: str = Field(..., description="The SMILES representation of the molecule.")

# --- Tool: Count Rotatable Bonds ---
@register_basic_tool(
    metadata=BasicToolMetadata(
        name="count_rotatable_bonds",
        description="Counts the number of rotatable bonds in a molecule, which influences its conformational flexibility.",
        tags=["rdkit", "chemistry", "descriptor", "structure"],
        author="RDKit Wrapper",
        parameters=CountRotatableBondsInput.model_json_schema()
    )
)
def count_rotatable_bonds(input: CountRotatableBondsInput) -> int:
    mol = _get_mol(input.smiles_string)
    return Descriptors.NumRotatableBonds(mol)

class CountHBondDonorsAcceptorsInput(BaseModel):
    smiles_string: str = Field(..., description="The SMILES representation of the molecule.")

# --- Tool: Count H-Bond Donors and Acceptors ---
@register_basic_tool(
    metadata=BasicToolMetadata(
        name="count_h_bond_donors_acceptors",
        description="Counts the number of hydrogen bond donors and acceptors in a molecule according to Lipinski's rules.",
        tags=["rdkit", "chemistry", "descriptor", "h-bond"],
        author="RDKit Wrapper",
        parameters=CountHBondDonorsAcceptorsInput.model_json_schema()
    )
)
def count_h_bond_donors_acceptors(input: CountHBondDonorsAcceptorsInput) -> dict:
    mol = _get_mol(input.smiles_string)
    donors = Lipinski.NumHDonors(mol)
    acceptors = Lipinski.NumHAcceptors(mol)
    return {"donors": donors, "acceptors": acceptors}

class CalculateMolecularSimilarityInput(BaseModel):
    smiles1: str = Field(..., description="SMILES string of the first molecule.")
    smiles2: str = Field(..., description="SMILES string of the second molecule.")
    metric: Literal["tanimoto", "dice", "cosine"] = Field("tanimoto", description="Similarity metric to use.")
    fingerprint_type: Literal["morgan", "rdkit"] = Field("morgan", description="Type of fingerprint to generate.")

# --- Tool: Calculate Molecular Similarity ---
@register_basic_tool(
    metadata=BasicToolMetadata(
        name="calculate_molecular_similarity",
        description="Calculates the molecular similarity between two molecules using specified fingerprint and metric.",
        tags=["rdkit", "chemistry", "similarity", "fingerprint"],
        author="RDKit Wrapper",
        parameters=CalculateMolecularSimilarityInput.model_json_schema()
    )
)
def calculate_molecular_similarity(input: CalculateMolecularSimilarityInput) -> float:
    mol1 = _get_mol(input.smiles1)
    mol2 = _get_mol(input.smiles2)

    if input.fingerprint_type == "morgan":
        fp1 = AllChem.GetMorganFingerprintAsBitVect(mol1, 2, nBits=2048)
        fp2 = AllChem.GetMorganFingerprintAsBitVect(mol2, 2, nBits=2048)
    elif input.fingerprint_type == "rdkit":
        fp1 = Chem.RDKFingerprint(mol1)
        fp2 = Chem.RDKFingerprint(mol2)
    else:
        raise ValueError(f"Unsupported fingerprint type: {input.fingerprint_type}. Choose 'morgan' or 'rdkit'.")

    if input.metric == "tanimoto":
        return DataStructs.TanimotoSimilarity(fp1, fp2)
    elif input.metric == "dice":
        return DataStructs.DiceSimilarity(fp1, fp2)
    elif input.metric == "cosine":
        return DataStructs.CosineSimilarity(fp1, fp2)
    else:
        raise ValueError(f"Unsupported similarity metric: {input.metric}. Choose 'tanimoto', 'dice', or 'cosine'.")

class FindSimilarMoleculesInListInput(BaseModel):
    query_smiles: str = Field(..., description="SMILES string of the query molecule.")
    smiles_list: List[str] = Field(..., description="A list of SMILES representations of molecules to search within.")
    threshold: float = Field(0.7, description="Similarity threshold (0.0 to 1.0).")
    metric: Literal["tanimoto", "dice", "cosine"] = Field("tanimoto", description="Similarity metric to use.")
    fingerprint_type: Literal["morgan", "rdkit"] = Field("morgan", description="Type of fingerprint to generate.")

# --- Tool: Find Similar Molecules in List ---
@register_basic_tool(
    metadata=BasicToolMetadata(
        name="find_similar_molecules_in_list",
        description="Finds molecules in a list that are similar to a query molecule based on a specified similarity threshold, fingerprint, and metric.",
        tags=["rdkit", "chemistry", "similarity", "screening"],
        author="RDKit Wrapper",
        parameters=FindSimilarMoleculesInListInput.model_json_schema()
    )
)
def find_similar_molecules_in_list(input: FindSimilarMoleculesInListInput) -> List[str]:
    if not (0.0 <= input.threshold <= 1.0):
        raise ValueError("Threshold must be between 0.0 and 1.0.")

    similar_molecules = []
    for target_smiles in input.smiles_list:
        try:
            similarity = calculate_molecular_similarity(
                CalculateMolecularSimilarityInput(
                    smiles1=input.query_smiles,
                    smiles2=target_smiles,
                    metric=input.metric,
                    fingerprint_type=input.fingerprint_type
                )
            )
            if similarity >= input.threshold:
                similar_molecules.append(target_smiles)
        except ValueError:
            # Skip invalid SMILES strings in the list or query
            continue
    return similar_molecules

class GenerateSimilarityMatrixReportInput(BaseModel):
    smiles_list: List[str] = Field(..., description="A list of SMILES representations of the molecules.")
    metric: Literal["tanimoto", "dice", "cosine"] = Field("tanimoto", description="Similarity metric to use.")
    fingerprint_type: Literal["morgan", "rdkit"] = Field("morgan", description="Type of fingerprint to generate.")
    output_format: Literal["csv", "json"] = Field("csv", description="The desired output format (csv or json).")

# --- Tool: Generate Similarity Matrix Report ---
@register_basic_tool(
    metadata=BasicToolMetadata(
        name="generate_similarity_matrix_report",
        description="Generates a similarity matrix report (CSV or JSON) for a list of molecules.",
        tags=["rdkit", "chemistry", "similarity", "report", "file-generation"],
        author="RDKit Wrapper",
        parameters=GenerateSimilarityMatrixReportInput.model_json_schema()
    )
)
def generate_similarity_matrix_report(input: GenerateSimilarityMatrixReportInput) -> str:
    if not input.smiles_list:
        return ""

    if input.output_format not in ["csv", "json"]:
        raise ValueError("Unsupported output format. Choose 'csv' or 'json'.")

    # Generate fingerprints for all valid molecules
    valid_mols = []
    valid_smiles = []
    for smiles in input.smiles_list:
        try:
            mol = _get_mol(smiles)
            valid_mols.append(mol)
            valid_smiles.append(smiles)
        except ValueError:
            # Skip invalid SMILES strings
            continue

    if not valid_mols:
        return ""

    fingerprints = []
    for mol in valid_mols:
        if input.fingerprint_type == "morgan":
            fingerprints.append(AllChem.GetMorganFingerprintAsBitVect(mol, 2, nBits=2048))
        elif input.fingerprint_type == "rdkit":
            fingerprints.append(Chem.RDKFingerprint(mol))
        else:
            raise ValueError(f"Unsupported fingerprint type: {input.fingerprint_type}. Choose 'morgan' or 'rdkit'.")

    # Calculate similarity matrix
    num_mols = len(valid_mols)
    similarity_matrix = [[0.0] * num_mols for _ in range(num_mols)]

    for i in range(num_mols):
        for j in range(i, num_mols):
            sim = 0.0
            if input.metric == "tanimoto":
                sim = DataStructs.TanimotoSimilarity(fingerprints[i], fingerprints[j])
            elif input.metric == "dice":
                sim = DataStructs.DiceSimilarity(fingerprints[i], fingerprints[j])
            elif input.metric == "cosine":
                sim = DataStructs.CosineSimilarity(fingerprints[i], fingerprints[j])
            else:
                raise ValueError(f"Unsupported similarity metric: {input.metric}. Choose 'tanimoto', 'dice', or 'cosine'.")
            similarity_matrix[i][j] = sim
            similarity_matrix[j][i] = sim # Matrix is symmetric

    if input.output_format == "csv":
        # Prepare CSV content
        header = ["SMILES"] + valid_smiles
        csv_lines = [",".join(header)]
        for i in range(num_mols):
            row = [valid_smiles[i]] + [f"{val:.4f}" for val in similarity_matrix[i]]
            csv_lines.append(",".join(row))
        return "\n".join(csv_lines)
    elif input.output_format == "json":
        # Prepare JSON content
        json_data = {
            "smiles": valid_smiles,
            "similarity_matrix": similarity_matrix
        }
        return json.dumps(json_data, indent=2)

class GenerateSimilarityHeatmapDataInput(BaseModel):
    smiles_list: List[str] = Field(..., description="A list of SMILES representations of the molecules.")
    metric: Literal["tanimoto", "dice", "cosine"] = Field("tanimoto", description="Similarity metric to use.")
    fingerprint_type: Literal["morgan", "rdkit"] = Field("morgan", description="Type of fingerprint to generate.")

# --- Tool: Generate Similarity Heatmap Data ---
@register_basic_tool(
    metadata=BasicToolMetadata(
        name="generate_similarity_heatmap_data",
        description="Generates data for a similarity heatmap (JSON) for a list of molecules.",
        tags=["rdkit", "chemistry", "similarity", "visualization", "heatmap"],
        author="RDKit Wrapper",
        parameters=GenerateSimilarityHeatmapDataInput.model_json_schema()
    )
)
def generate_similarity_heatmap_data(input: GenerateSimilarityHeatmapDataInput) -> dict:
    if not input.smiles_list:
        return {"labels": [], "matrix": [], "min_val": 0.0, "max_val": 0.0}

    # Generate fingerprints for all valid molecules
    valid_mols = []
    valid_smiles = []
    for smiles in input.smiles_list:
        try:
            mol = _get_mol(smiles)
            valid_mols.append(mol)
            valid_smiles.append(smiles)
        except ValueError:
            # Skip invalid SMILES strings
            continue

    if not valid_mols:
        return {"labels": [], "matrix": [], "min_val": 0.0, "max_val": 0.0}

    fingerprints = []
    for mol in valid_mols:
        if input.fingerprint_type == "morgan":
            fingerprints.append(AllChem.GetMorganFingerprintAsBitVect(mol, 2, nBits=2048))
        elif input.fingerprint_type == "rdkit":
            fingerprints.append(Chem.RDKFingerprint(mol))
        else:
            raise ValueError(f"Unsupported fingerprint type: {input.fingerprint_type}. Choose 'morgan' or 'rdkit'.")

    # Calculate similarity matrix
    num_mols = len(valid_mols)
    similarity_matrix = [[0.0] * num_mols for _ in range(num_mols)]
    min_val = 1.0
    max_val = 0.0

    for i in range(num_mols):
        for j in range(i, num_mols):
            sim = 0.0
            if input.metric == "tanimoto":
                sim = DataStructs.TanimotoSimilarity(fingerprints[i], fingerprints[j])
            elif input.metric == "dice":
                sim = DataStructs.DiceSimilarity(fingerprints[i], fingerprints[j])
            elif input.metric == "cosine":
                sim = DataStructs.CosineSimilarity(fingerprints[i], fingerprints[j])
            else:
                raise ValueError(f"Unsupported similarity metric: {input.metric}. Choose 'tanimoto', 'dice', or 'cosine'.")
            similarity_matrix[i][j] = sim
            similarity_matrix[j][i] = sim # Matrix is symmetric

            if sim < min_val:
                min_val = sim
            if sim > max_val:
                max_val = sim

    return {
        "labels": valid_smiles,
        "matrix": similarity_matrix,
        "min_val": min_val,
        "max_val": max_val
    }
