import joblib
import pandas as pd
import numpy as np
from Bio import Entrez, Medline
from sklearn.feature_extraction.text import TfidfVectorizer
from sklearn.svm import SVC
from imblearn.pipeline import Pipeline
from imblearn.over_sampling import SMOTE
from sklearn.model_selection import train_test_split
from sklearn.metrics import accuracy_score, classification_report, confusion_matrix
import logging
import os

# Configure logging
logging.basicConfig(filename="drug_model_error.log", level=logging.ERROR)

# Constants (for clarity and easier modification)
MODEL_PATH = "interaction_classifier.pkl"
VECTORIZER_PATH = "vectorizer.pkl"
EMAIL = "sannody0@outlook.sa"  # your email
DEFAULT_SEARCH_TERM = "drug interactions"
DEFAULT_RETMAX = 1000

def fetch_pubmed_articles(term, retmax=DEFAULT_RETMAX):
    """Fetch PubMed articles with robust error handling."""
    Entrez.email = EMAIL
    try:
        handle = Entrez.esearch(db="pubmed", term=term, retmax=retmax)
        record = Entrez.read(handle)
        id_list = record["IdList"]
        handle = Entrez.efetch(db="pubmed", id=id_list, rettype="medline", retmode="text")
        records = Medline.parse(handle)
        return list(records)
    except Exception as e:
        logging.error(f"Error fetching PubMed articles: {e}")
        raise  # Re-raise the exception after logging

def prepare_data_from_pubmed(term, retmax=DEFAULT_RETMAX):
    """Prepare training data from PubMed articles."""
    articles = fetch_pubmed_articles(term, retmax)
    data = [{"abstract": record.get("AB", ""), "interaction": int("interaction" in record.get("TI", "").lower())} for record in articles]
    df = pd.DataFrame(data)
    df.dropna(subset=["abstract"], inplace=True)

    vectorizer = TfidfVectorizer(stop_words="english")
    X = vectorizer.fit_transform(df["abstract"])

    return X, df['interaction'], vectorizer  # Return labels as well

def train_and_save_model(X, y, model_filename=MODEL_PATH, vectorizer_filename=VECTORIZER_PATH):
    """Train, evaluate, and save the model."""
    X_train, X_test, y_train, y_test = train_test_split(X, y, test_size=0.2, random_state=42)

    # Check if there's enough data for training
    if len(np.unique(y_train)) <= 1:
        logging.error("Not enough data for training. Need at least two classes.")
        raise ValueError("Insufficient training data")

    model = Pipeline([
        ("sampling", SMOTE(random_state=42)),
        ("classifier", SVC(probability=True))
    ])

    model.fit(X_train, y_train)
    y_pred = model.predict(X_test)

    print("Accuracy:", accuracy_score(y_test, y_pred))
    print("Classification Report:\n", classification_report(y_test, y_pred))
    print("Confusion Matrix:\n", confusion_matrix(y_test, y_pred))

    joblib.dump(model, model_filename)
    joblib.dump(vectorizer, vectorizer_filename)
    return model, vectorizer

def load_model(model_path=MODEL_PATH, vectorizer_path=VECTORIZER_PATH):
    """Load the trained model and vectorizer."""
    model = joblib.load(model_path)
    vectorizer = joblib.load(vectorizer_path)
    return model, vectorizer

def predict_interactions(model, vectorizer, articles):
    df = pd.DataFrame(articles)
    X = vectorizer.transform(df['abstract'])

    # طباعة عدد الميزات للتحقق من التناسق
    print(f"Number of features in vectorizer: {X.shape[1]}")
    print(f"Expected number of features: {vectorizer.get_feature_names_out().shape[0]}")

    try:
        # Try to use predict_proba if available
        y_prob = model.predict_proba(X)[:, 1]
    except AttributeError as e:
        logging.warning(f"Model does not support predict_proba: {e}")
        y_prob = [None] * X.shape[0]  # Assign None or 0.5 as default probability

    predictions = model.predict(X)
    return predictions, y_prob

# Main function for training (if needed)
if __name__ == "__main__":
    if not os.path.exists(MODEL_PATH) or not os.path.exists(VECTORIZER_PATH):
        print("Pre-trained model not found. Training a new model...")
        X, y, vectorizer = prepare_data_from_pubmed(DEFAULT_SEARCH_TERM, retmax=DEFAULT_RETMAX)
        train_and_save_model(X, y)
    else:
        print("Loaded pre-trained model and vectorizer.")
