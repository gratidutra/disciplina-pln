import os
import re
from tqdm import tqdm
import spacy
import pandas as pd

# Carregue sua lista NCBI gene
def load_ncbi_gene_names(filepath):
    gene_names = set()
    with open(filepath, "r", encoding="utf-8") as f:
        for line in f:
            if line.startswith("#"):
                continue
            cols = line.strip().split("\t")
            if len(cols) < 3:
                continue
            symbol = cols[2].strip()
            if symbol:
                gene_names.add(symbol.upper())
    return gene_names

# Heurística para eliminar falsos positivos
def is_valid_gene(ent):
    if re.search(r"[,]|M$", ent):
        return False
    if "[" in ent or "]" in ent:
        return False
    if "." in ent and not re.fullmatch(r"[A-Z]\.", ent):
        return False
    if "-" in ent or " " in ent:
        return False
    if len(ent) <= 2:
        return False
    return True

BLACKLIST = {
    "e. coli", "gdn", "gdn hcl", "gndhcl", "28,880 m", "bldg. 3",
    "solvent-membrane", "trp", "tyr", "omp", "omps", "s2", "pagp",
    "vdac", "foma", "yop", "transmembrane α -helices"
}

# Verifica se gene está na lista NCBI
def is_known_gene(ent):
    return ent.upper() in valid_genes

# Função de extração — seu código de extração NER com scispaCy e entity ruler
def extract_entities(model, document_text, document_name):
    nlp = model.load()

    if "entity_ruler" in nlp.pipe_names:
        ruler = nlp.get_pipe("entity_ruler")
        ruler.clear()
    else:
        ruler = nlp.add_pipe("entity_ruler", after="ner")

    patterns = [
        # Padrões MUTATION
        {"label": "MUTATION", "pattern": [{"TEXT": {"REGEX": r"^(Δ|[A-Z])[0-9]+(?:[A-Za-z]*|del|ins|dup|fs|X)$"}}]},
        {"label": "MUTATION", "pattern": [{"LOWER": "mutation"}, {"TEXT": {"REGEX": r"^(Δ|[A-Z])[0-9]+(?:[A-Za-z]*|del|ins|dup|fs|X)$"}}]},
        {"label": "MUTATION", "pattern": [{"TEXT": {"REGEX": r"^p\.[A-Z][a-z]{2}[0-9]+(?:[A-Za-z]*|del|ins|dup)"}}]},
        {"label": "MUTATION", "pattern": [{"TEXT": {"REGEX": r"^c\.[0-9]+(?:_[0-9]+)?(?:del|ins|dup)[A-Z]*"}}]},

        # Padrões METHOD
        {"label": "METHOD", "pattern": [{"LOWER": "urea"}]},
        {"label": "METHOD", "pattern": [{"LOWER": "gdn"}, {"LOWER": "hcl"}]},
        {"label": "METHOD", "pattern": [{"LOWER": "guanidine"}, {"LOWER": "hydrochloride"}]},
        {"label": "METHOD", "pattern": [{"LOWER": "thermal"}, {"LOWER": "shift"}, {"LOWER": "assay"}]},
        {"label": "METHOD", "pattern": [{"LOWER": "tsa"}]},

        # Padrões MEASURE
        {"label": "MEASURE", "pattern": [{"LOWER": "differential"}, {"LOWER": "scanning"}, {"LOWER": "calorimetry"}]},
        {"label": "MEASURE", "pattern": [{"LOWER": "dsc"}]},
        {"label": "MEASURE", "pattern": [{"LOWER": "circular"}, {"LOWER": "dichroism"}, {"LOWER": "spectroscopy"}]},
        {"label": "MEASURE", "pattern": [{"LOWER": "cd"}, {"LOWER": "spectroscopy"}]},
        {"label": "MEASURE", "pattern": [{"LOWER": "cd"}]}
    ]

    ruler.add_patterns(patterns)

    doc = nlp(document_text)

    custom_labels = {"MUTATION", "METHOD", "MEASURE", "GENE_OR_GENE_PRODUCT"}
    seen = set()
    entity = []
    label = []

    for ent in doc.ents:
        if ent.label_ in custom_labels and ent.text.lower() not in seen:
            entity.append(ent.text)
            label.append(ent.label_)
            seen.add(ent.text.lower())

    return {
        "pmid": document_name,
        "entity": entity,
        "label": label
    }

class SpacyModel:
    def load(self):
        return spacy.load("en_ner_bionlp13cg_md")

# --- EXECUÇÃO PRINCIPAL ---

# Carrega lista NCBI Gene (antes do loop)
ncbi_gene_path = "data/Homo_sapiens.gene_info"  # Ajuste o caminho conforme seu arquivo
valid_genes = load_ncbi_gene_names(ncbi_gene_path)

caminho_pasta = "data/md_cleaned"
model = SpacyModel()
todos_resultados = []
arquivos_md = [f for f in os.listdir(caminho_pasta) if f.endswith(".md")]

for nome_arquivo in tqdm(arquivos_md, desc="Processando arquivos"):
    caminho_completo = os.path.join(caminho_pasta, nome_arquivo)
    with open(caminho_completo, "r", encoding="utf-8") as f:
        document_text = f.read()

    resultado = extract_entities(model, document_text, nome_arquivo)

    # Aqui aplica o filtro *somente* para genes, aceita os outros labels direto
    for ent, label in zip(resultado["entity"], resultado["label"]):
        ent_clean = ent.strip().lower()
        if label == "GENE_OR_GENE_PRODUCT":
            if ent_clean in BLACKLIST or not is_valid_gene(ent) or not is_known_gene(ent):
                continue  # Ignora falso positivo
        # Salva MUTATION, METHOD, MEASURE direto, sem filtro extra
        todos_resultados.append({
            "pmid": resultado["pmid"],
            "entity": ent,
            "label": label,
        })

df = pd.DataFrame(todos_resultados)

#df.to_csv("result_model.csv", index=False)

# Função que junta as entidades por label para cada pmid
df_wide = df.groupby(['pmid', 'label'])['entity'] \
    .apply(lambda x: ', '.join(sorted(set(x)))) \
    .unstack(fill_value='')

# Se quiser resetar índice para virar coluna:
df_wide = df_wide.reset_index()

# Salvar o novo CSV
df_wide.to_csv('result_model.csv', index=False)
