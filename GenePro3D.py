import streamlit as st
from stmol import showmol
import py3Dmol
import requests
import biotite.structure.io as bsio

# Function to transcribe a DNA sequence into an mRNA sequence
def transcribe(dna_sequence):
    return dna_sequence.replace('T', 'U')

# Function to translate an mRNA sequence into an amino acid sequence
def translate(mrna_sequence):
    codon_table = {
        'AUG': 'M', 'UUU': 'F', 'UUC': 'F', 'UUA': 'L', 'UUG': 'L', 'CUU': 'L', 'CUC': 'L', 'CUA': 'L', 'CUG': 'L',
        'AUU': 'I', 'AUC': 'I', 'AUA': 'I', 'GUU': 'V', 'GUC': 'V', 'GUA': 'V', 'GUG': 'V',
        'UCU': 'S', 'UCC': 'S', 'UCA': 'S', 'UCG': 'S', 'CCU': 'P', 'CCC': 'P', 'CCA': 'P', 'CCG': 'P',
        'ACU': 'T', 'ACC': 'T', 'ACA': 'T', 'ACG': 'T', 'GCU': 'A', 'GCC': 'A', 'GCA': 'A', 'GCG': 'A',
        'UAU': 'Y', 'UAC': 'Y', 'CAU': 'H', 'CAC': 'H', 'CAA': 'Q', 'CAG': 'Q', 'AAU': 'N', 'AAC': 'N',
        'AAA': 'K', 'AAG': 'K', 'GAU': 'D', 'GAC': 'D', 'GAA': 'E', 'GAG': 'E', 'UGU': 'C', 'UGC': 'C',
        'UGG': 'W', 'CGU': 'R', 'CGC': 'R', 'CGA': 'R', 'CGG': 'R', 'AGU': 'S', 'AGC': 'S', 'AGA': 'R', 'AGG': 'R',
        'UAA': '*', 'UAG': '*', 'UGA': '*'
    }
    protein = ''
    start_codon = 'AUG'
    stop_codons = ['UAA', 'UAG', 'UGA']
    start_index = mrna_sequence.find(start_codon)
    if start_index == -1:
        return protein

    i = start_index
    while i + 3 <= len(mrna_sequence):
        codon = mrna_sequence[i:i + 3]
        if codon in stop_codons:
            break
        amino_acid = codon_table.get(codon, 'X')  # 'X' for an unknown codon
        protein += amino_acid
        i += 3

    return protein
    pass

# Function to render protein 3D structure
def render_mol(pdb):
    pdbview = py3Dmol.view()
    pdbview.addModel(pdb, 'pdb')
    pdbview.setStyle({'cartoon': {'color': 'spectrum'}})
    pdbview.setBackgroundColor('white')
    pdbview.zoomTo()
    pdbview.zoom(2, 800)
    pdbview.spin(True)
    showmol(pdbview, height = 500,width=800)

# Streamlit web application
# Function to predict 3D structure
def update(sequence):
    headers = {
        'Content-Type': 'application/x-www-form-urlencoded',
    }
    response = requests.post('https://api.esmatlas.com/foldSequence/v1/pdb/', headers=headers, data=sequence, timeout=60)
    name = sequence[:3] + sequence[-3:]
    pdb_string = response.content.decode('utf-8')

    with open('predicted.pdb', 'w') as f:
        f.write(pdb_string)
    try:
    	struct = bsio.load_structure('predicted.pdb', extra_fields=["b_factor"])
    	b_value = round(struct.b_factor.mean(), 4)

    	# Display protein structure
    	st.subheader('Predicted Protein 3D Structure')
    	render_mol(pdb_string)

    	# plDDT value is stored in the B-factor field
    	st.subheader('plDDT')
    	st.write('plDDT is a per-residue estimate of the confidence in prediction on a scale from 0-100.')
    	st.info(f'plDDT: {b_value}')

    	st.download_button(
        	label="Download PDB",
        	data=pdb_string,
        	file_name='predicted.pdb',
        	mime='text/plain',
    	)
    except ValueError as e:
        # Display custom error message when no models are found
        st.error("ValueError: The file has 0 models, the given model 1 does not exist")
    
# Define the render_mol function for rendering protein 3D structure
def render_mol(pdb):
    pdbview = py3Dmol.view()
    pdbview.addModel(pdb, 'pdb')
    pdbview.setStyle({'cartoon': {'color': 'spectrum'}})
    pdbview.setBackgroundColor('white')
    pdbview.zoomTo()
    pdbview.zoom(2, 800)
    pdbview.spin(True)
    showmol(pdbview, height = 500,width=800)

# Sidebar with information
st.set_page_config(layout = 'wide')
with st.sidebar:
    st.subheader("GenePro3D")
    st.markdown("This is a sequence translator and 3D protein structure predictor. It transcribes DNA sequences and translates the resulting mRNA sequences into protein sequences before predicting the 3D structure of the protein. Its protein 3D structure prediction is based on the [ESM](https://esmatlas.com/about)-2 language model. For additional information, read [news article](https://www.nature.com/articles/d41586-022-03539-1) and [research article](https://www.biorxiv.org/content/10.1101/2022.07.20.500902v2).")

# Add 10 blank lines using HTML <br>
    for _ in range(10):
        st.write("<br>", unsafe_allow_html=True)
    
    st.subheader("Creators")
    st.markdown("[Okello Harrison Onyango](https://www.researchgate.net/profile/Harrison-Onyango)  \n[Masinde Muliro University of Science and Technology](https://mmust.ac.ke/)")

st.title("GenePro3D Predictor")

# User input for DNA or Protein sequence
st.markdown("**<p style='font-size: 20px; font-weight: bold;'>Select the Type of Sequence:</p>**", unsafe_allow_html=True)
sequence_type = st.radio("", ("DNA", "Protein"))


if sequence_type == "DNA":
    # User input for DNA sequence
    st.markdown("**Enter a DNA Sequence (e.g., ATGCCGAT...):**")
    dna_sequence = st.text_area("", height=100).replace(" ", "").upper()

    if st.button("Translate"):
        # Transcribe the DNA sequence into mRNA
        mrna_sequence = transcribe(dna_sequence)

        # Translate the mRNA sequence into an amino acid sequence
        amino_acid_sequence = translate(mrna_sequence)
        st.markdown("**Protein Sequence:**")
        st.write(amino_acid_sequence)

        # Predict 3D structure and render it
        update(sequence=amino_acid_sequence)
        pass

elif sequence_type == "Protein":
    # User input for Protein sequence
    st.markdown("**Enter a Protein Sequence (e.g., MGSSHHHHHH...):**")
    protein_sequence = st.text_area("", height=100)

    if st.button("Predict 3D Structure"):
        # Predict 3D structure and render it
        update(sequence=protein_sequence)
else:
    st.warning("Select a Sequence Type.")


