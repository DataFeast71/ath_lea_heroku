import streamlit as st
import altair as alt
import pandas as pd
import json
from collections import Counter

st.title("Arabidopsis LEA Display information Web App")

st.markdown("""
This app retrieves the information from *Arabidopsis thaliana* LEA (late embriogenesis abundant) proteins and the results from different predictors:
* VSL2
* IUPRED
* ANCHOR
* MOBIDB
* FUZZPRED

All the LEAS protein from *A. thaliana* are present in the JSON file.
""")

st.sidebar.header("Use Input Features")

# Load Data from JSON file


def load_data(file):
    with open(file, "r") as f:
        data = json.load(f)

    # Data preparation
    uniprot = [protein["uniprot"] for _, protein in data.items()]
    tair_protein = [protein["tair_protein"] for _, protein in data.items()]
    tair_gen = [protein["tair_gen"] for _, protein in data.items()]
    classPFAM = [protein["classification"]["ClassPFMA"]
                 for _, protein in data.items()]
    classAcovs = [protein["classification"]["ClassAcovs"]
                  for _, protein in data.items()]
    pfam = [protein["domains"]["PFAM"] for _, protein in data.items()]
    gravy = [protein["parameters"]["GRAVY"] for _, protein in data.items()]
    length = [len(protein["sequence"]) for _, protein in data.items()]

    # Make a Dataframe
    df = pd.DataFrame({"UNIPROT": uniprot,
                       "TAIR-Protein": tair_protein,
                       "TAIR-Gen": tair_gen,
                       "Class-PFAM": classPFAM,
                       "Class-Acovs": classAcovs,
                       "PFAM": pfam,
                       "GRAVY": gravy,
                       "Length": length})
    return data, df


data, df_leas = load_data("./AthLeasDB.json")

# Side Bar - Class Selection
class_selected = list(df_leas['Class-Acovs'].unique())
selected_class = st.sidebar.multiselect(
    "Groups by Acovs", class_selected, class_selected)

# Filtering Data by Class
df_selected = df_leas[df_leas['Class-Acovs'].isin(selected_class)]

# Side Bar - Protein selection
protein_selected = list(df_selected['TAIR-Protein'].unique())
selected_protein = st.sidebar.selectbox("Protein selection", protein_selected)

# Side Bar- Predictor selection
predictors_selected = ["FUZZPRED", "VSL2", "ANCHOR", "IUPRED", "MOBIDB"]
selected_predictor = st.sidebar.selectbox(
    "Predictor to plot:", predictors_selected)

st.header("Display LEA information from the Group(s) selected")
st.write("Data Dimension: {} rows and {} columns.".format(
    df_selected.shape[0], df_selected.shape[1]))

st.dataframe(df_selected)

#  Download Data? a feature?

#####################
# Counting Aminoacids
#####################
# Counting Aminoacids residues for one protein
st.subheader("Aminoacid Residue Count.")

AA_count = Counter(data[selected_protein]['sequence'])

# Display in a DataFrame
df_aa = pd.DataFrame.from_dict(AA_count,
                               orient="index").rename({
                                   0: "Count"}, axis=1).reset_index()
df_aa = df_aa.rename({'index': "Aminoacid"}, axis=1)
df_aa['Percentage'] = df_aa['Count'].apply(lambda x: 100 * x/df_aa.Count.sum())
df_aa = df_aa.sort_values("Count", ascending=False)

st.dataframe(df_aa)

# Display Altair chart
st.subheader("Barplot Aminiacid Counting.")
p = alt.Chart(df_aa).mark_bar().encode(
    x="Count",
    y=alt.Y("Aminoacid", sort="-x")
)
st.altair_chart(p)

########################
# Predictor results
########################
# Area Chart for predictor
df_predictor = pd.DataFrame(
    {"Scores": data[selected_protein]["predictors"][selected_predictor]['scores']})
df_predictor['Pos'] = list(range(1, df_predictor.shape[0] + 1))

st.subheader("Results form {} predictor".format(selected_predictor))
p = alt.Chart(df_predictor).mark_area(
    line={"color": 'darkblue'},
    color=alt.Gradient(
        gradient="linear",
        stops=[alt.GradientStop(color='white', offset=0),
               alt.GradientStop(color='darkblue', offset=1)],
        x1=1,
        x2=1,
        y1=1,
        y2=0
    )
).encode(
    alt.X("Pos", title="Position"),
    alt.Y("Scores")
).properties(
    width=600,
    height=180
)
st.altair_chart(p)
