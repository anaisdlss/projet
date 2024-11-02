import streamlit as st
from Bio.Seq import Seq
from Bio import SeqIO
from io import StringIO
import re

# Définition des catégories d'acides aminés avec leurs couleurs
hydrophobic_amino_acids = "GAILMFVWP"
polar_amino_acids = "STCNQY"
positive_amino_acids = "KRH"
negative_amino_acids = "DE"


def color_protein_sequence(protein_sequence):
    """Colorie les acides aminés dans une séquence protéique.

    Chaque acide aminé est coloré en fonction de sa catégorie :
    - Hydrophobe : jaune
    - Polaire : vert
    - Positif : rouge
    - Négatif : bleu

    Parameters
    ----------
    protein_sequence : str
        La séquence protéique à colorer.

    Returns
    -------
    str
        La séquence colorée au format HTML.
    """
    colored_sequence = ""
    for aa in protein_sequence:
        if aa in hydrophobic_amino_acids:
            colored_sequence += f'<span style="color:yellow">{aa}</span>'
        elif aa in polar_amino_acids:
            colored_sequence += f'<span style="color:green">{aa}</span>'
        elif aa in positive_amino_acids:
            colored_sequence += f'<span style="color:red">{aa}</span>'
        elif aa in negative_amino_acids:
            colored_sequence += f'<span style="color:blue">{aa}</span>'
        else:
            colored_sequence += aa

    return colored_sequence


def translate_dna(sequence, frame):
    """Traduit une séquence d'ADN en protéine.

    Cette fonction prend une séquence d'ADN et un cadre de lecture,
    et retourne la séquence de protéine traduite.

    Parameters
    ----------
    sequence : str
        La séquence d'ADN à traduire.
    frame : int
        Le cadre de lecture (1, 2, ou 3).

    Returns
    -------
    str or None
        La séquence de protéine traduite ou None en cas d'erreur.
    """
    try:
        # Ajustement du cadre de lecture (Streamlit utilise 1, 2, 3, mais BioPython commence à 0)
        dna_seq = Seq(sequence[frame - 1:])
        protein_seq = dna_seq.translate(to_stop=False)
        return str(protein_seq)
    except Exception as e:
        st.error(f"Erreur lors de la traduction : {e}")
        return None


def readfasta(fasta):
    """Lit un fichier FASTA et retourne une liste de séquences.

    La fonction traite le contenu d'un fichier FASTA pour extraire
    les noms de séquences et les séquences correspondantes.

    Parameters
    ----------
    fasta : str
        Le contenu d'un fichier FASTA sous forme de chaîne.

    Returns
    -------
    list of tuple
        Une liste de tuples contenant les noms de séquences (str) et
        les séquences correspondantes (str).
    """
    lignes = fasta.splitlines()
    seq = ""
    L = []
    seqnom = ""
    for ligne in lignes:
        ligne = ligne.strip()
        if ligne.startswith(">"):  # Si c'est un nom de séquence
            if seqnom != "":
                L.append((seqnom, seq))
            seqnom = ligne  # On met à jour le nom de la nouvelle séquence
            seq = ""  # On réinitialise la séquence
        else:
            seq += ligne  # On ajoute la ligne à la séquence
    # Ajouter la dernière séquence après la boucle
    if seqnom != "":
        L.append((seqnom, seq))
    return L


# Interface de l'application
st.title("Traduction de séquences d'ADN en protéines")

# Formulaire d'entrée de séquences ADN au format FASTA
st.header("Entrer des séquences nucléiques au format FASTA")
fasta_input = st.text_area(
    "Saisissez vos séquences nucléiques au format FASTA ici")

# Formulaire pour envoyer un fichier
fichiermulti = st.file_uploader(
    "Choisissez un fichier fasta", accept_multiple_files=True)
Ltextfichier = []
for fichier in fichiermulti:
    textfichier = fichier.read().decode("utf-8")
    Ltextfichier.append(textfichier)
    # Forme une longue chaîne de caractères comme si on avait input directement dans la zone texte
    finalfichier = "\n".join(Ltextfichier)

# Choix du cadre de lecture
frame = st.selectbox(
    "Choisissez le cadre de lecture (1, 2, ou 3)", options=[1, 2, 3])

# Bouton pour effectuer la traduction
translated_sequences = []
Lsequence = []
if st.button("Traduire"):
    if fasta_input:
        # Lecture des séquences FASTA fournies par l'utilisateur
        Lsequence = readfasta(fasta_input)

    if fichiermulti:
        Lsequence = readfasta(finalfichier)
    else:
        st.warning("Veuillez entrer des séquences nucléiques au format FASTA.")
    if len(Lsequence) < 0:
        for dna_seq in Lsequence:
            seq_name = dna_seq[0]
            protein_seq = translate_dna(dna_seq[1], frame)
            colored_protein = color_protein_sequence(protein_seq)
            translated_sequences.append(f"> {seq_name}\n{colored_protein}")


# Affichage des séquences traduites
st.header("Séquences protéiques traduites")
for seq in translated_sequences:
    st.markdown(seq, unsafe_allow_html=True)

# Générer le fichier FASTA pour le téléchargement
fasta_output = "\n".join([re.sub("<[^<]+?>", "", seq)
                         for seq in translated_sequences])
fasta_file = StringIO(fasta_output)

# Bouton pour télécharger le fichier FASTA
st.download_button(
    label="Télécharger les séquences traduites",
    data=fasta_file.getvalue(),
    file_name="translated_sequences.fasta",
    mime="text/plain",
)
